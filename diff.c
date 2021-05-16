//
// TODO 
//   better memory handling for library
//   split into library interface
//
#include <err.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef MIN
# define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

// Pair of serial (line number) and line hash
//
typedef struct serhash_t serhash_t;
struct serhash_t {
    unsigned serial;
    unsigned hash;
};

// Pair of serial (line number) and 'last' flag, indicating this
// is the last line in an equivalence class
//
typedef struct eclass_t eclass_t;
struct eclass_t {
    unsigned serial;
    int last;
};

// A generic vector 
//
typedef struct vec_t vec_t;
struct vec_t {
    size_t elesize;         // size of each element
    size_t alloc;           // # of elements allocated
    size_t n;               // # of elements in use
    char *data;             // -> storage
};

// A k-candidate, with line number for the left and right file each
// prev is used to link chains of k-candidates
//
typedef struct candidate_t candidate_t;
struct candidate_t {
    int a;
    int b;
    candidate_t *prev;
};

// A range of lines, [start..end)
//
typedef struct range_t range_t; 
struct range_t {
    int start;
    int end;
};

// A change between the two files; two ranges of lines that differ
//
typedef struct delta_t delta_t;
struct delta_t {
    range_t left;
    range_t right;
};

// Overall state of the algorithm with everything needed to clean up
// memory
//
typedef struct diffstate_t diffstate_t;
struct diffstate_t {
    FILE *fp[2];            // left and right file pointers 
    vec_t V;                // (hash, serial) pairs of right hand file
    vec_t E;                // (serial, last) pairs of right hand file eq classes
    vec_t P;                // left hand file points by line into matching eq class in E
    candidate_t **K;        // rightmost k-candidate pointers for each K
    int *J;                 // left-right matchines at end of alg
    vec_t dv;               // vector of delta ranges for output construction
};

/*****************************************************************************
 * Generic vector  
 */

//
// Initialize an empty vector
//
// Returns 0 on success or -ENOMEM on out of memory
//
int vec_init(vec_t *pv, size_t elesize)
{
    const unsigned INITIAL_ALLOC = 16;

    pv->elesize = elesize;
    pv->alloc = INITIAL_ALLOC;
    pv->n = 0;
    pv->data = malloc(elesize * INITIAL_ALLOC);
    if (pv->data == NULL) {
        return -ENOMEM;
    }
    return 0;
}

//
// Add an element to a vector
//
// Returns the new element or NULL on out of memory
//
void *vec_append(vec_t *pv)
{
    size_t oldsize, newsize;
    size_t newalloc;
    char *newp;

    if (pv->n == pv->alloc) {
        oldsize = pv->alloc * pv->elesize;
        newalloc *= 2;
        newsize = newalloc * pv->elesize;
        if (newsize <= oldsize) {
            // size has wrapped around, too big!
            return NULL;
        }
        newp = realloc(pv->data, newsize);
        if (newp == NULL) {
            newalloc = pv->alloc + 1;
            newsize = newalloc * pv->elesize;
            newp = realloc(pv->data, newsize);
            if (newp == NULL) {
                return NULL;
            }
        }
        pv->alloc = newalloc;
        pv->data = newp;
    }

    pv->n++;
    return &pv->data[pv->elesize * (pv->n - 1)];
}

// clamp a vector so its allocated size is the exactly the number of elements
// in use. used to recover space reserved for new elements when it's known
// that the vector will never again need to grow.
//
// Return 0 on success or -ENOMEM on out of memory
//
int vec_clamp(vec_t *pv)
{
    char *newp;

    newp = realloc(pv->data, pv->elesize * pv->n);
    if (newp == NULL) {
        // shouldn't happen since we're shrinking
        return -ENOMEM;
    }
    pv->alloc = pv->n;
    return 0;
}

//
// Free a vector's data
//
void vec_free(vec_t *v)
{
    if (v->data) {
        free(v->data);
        v->data = NULL;
        v->alloc = 0;
        v->n = 0;
    }
}

/*****************************************************************************
 * serial-hash vector
 */
#define SH_INIT(sh) vec_init(sh, sizeof(serhash_t))
#define SH_APPEND(sh) ((serhash_t*)vec_append(sh))
#define SH_AT(sh, idx) \
    ( ((serhash_t*)((sh)->data)) + (idx) )

/*****************************************************************************
 * equivalence class vector
 */
#define EV_INIT(ev) vec_init(ev, sizeof(eclass_t))
#define EV_APPEND(ev) ((eclass_t*)vec_append(ev))
#define EV_AT(ev, idx) \
    ( ((eclass_t*)((ev)->data)) + (idx) )

/*****************************************************************************
 * vector of integers
 */
#define INTV_INIT(iv) vec_init(iv, sizeof(int))
#define INTV_APPEND(iv) ((int*)(vec_append(iv)))
#define INTV_AT(iv, idx) \
    ( ((int*)((iv)->data)) + (idx) )

/*****************************************************************************
 * vector of deltas
 */
#define DV_INIT(dv) vec_init(dv, sizeof(delta_t))
#define DV_APPEND(dv) ((delta_t*)(vec_append(dv)))
#define DV_AT(dv, idx) \
    ( ((delta_t*)((dv)->data)) + (idx) )


/*****************************************************************************
 * serial-hash functions
 */

//
// Compare two serial-hash nodes on hash, then serial
//
int sh_compare(const void *l, const void *r)
{
    const serhash_t *ln = (const serhash_t *)l;
    const serhash_t *rn = (const serhash_t *)r;

    if (ln->hash < rn->hash) {
        return -1;
    }

    if (ln->hash > rn->hash) {
        return 1;
    }

    if (ln->serial < rn->serial) {
        return -1;
    }

    if (ln->serial > rn->serial) {
        return 1;
    }

    return 0;
}

//
// Sort a serial-hash vector on hash, then serial
//
void sh_sort(vec_t *sh)
{
    // NB we skip the first placeholder element
    qsort(sh->data + sizeof(serhash_t), sh->n - 1, sizeof(serhash_t), &sh_compare);
}

//
// Search for the last occurance of an equivalence class (based on hash)
// in the given sh vector. Return the index or 0 if there is no match.
//
unsigned sh_search(serhash_t *V, eclass_t *E, int n, unsigned hash)
{
    unsigned low = 1;
    unsigned high = n - 1;

    // since the array is sorted, we can binary search to find any class
    // member and then back up to find the first class element
    while (low <= high) {
        unsigned mid = (low + high) / 2;
        unsigned mid_h = V[mid].hash;

        if (hash < mid_h) {
            high = mid - 1;
        } else if (hash > mid_h) {
            low = mid + 1;
        } else {
            while (!E[mid-1].last) {
                mid--;
            }

            return mid;
        }
    }

    return 0;
}

/*****************************************************************************
 * candidates
 */

// constructor function for a k-candidate
//
candidate_t* candidate(int a, int b, candidate_t *prev)
{
    candidate_t *p = (candidate_t *)malloc(sizeof(candidate_t));
    if (p != NULL) {
        p->a = a;
        p->b = b;
        p->prev = prev;
    }
    return p;
}

/*****************************************************************************
 * file parsing and algorithm
 */

// If ch is an end of line char, skip any possible newline combination
// Also check for EOF due to error.
//
// Return 0 on success, -EIO on I/O error
//
int skipeol(FILE *fp, int ch) 
{
    if (ch == EOF && ferror(fp)) {
        return -EIO;
    }

    if (ch == '\r') {
        ch = fgetc(fp);
        if (ch != '\n') {
            if (ungetc(ch, fp) == EOF) {
                return -EIO;
            }
        } else if (ch == EOF && ferror(fp)) {
            return -EIO;
        }
    }
    return 0;
}

//
// Read the next line from 'fp' and compute its hash.
// 
// Return 1 on success with 'phash' containing the hash.
// Return 0 on EOF w/ no data to hash. 
// Return -EIO on I/O error.
//
int next_line(FILE *fp, unsigned *phash)
{
    int rc;
    int ch;
    int len = 0;
    const unsigned HASH0 = 5381;
    unsigned hash = HASH0;

    if (feof(fp)) {
        return 0;
    }

    while ((ch = fgetc(fp)) != EOF && ch != '\r' && ch != '\n') {
        len++;
        hash = ((hash << 5) + hash) + ch;
    }

    if ((rc = skipeol(fp, ch)) != 0) {
        return rc;
    }

    if (ch == EOF && len == 0) {
        // empty line at EOF
        //
        return 0;
    }

    *phash = hash;
    return 1;
}

// Skip the next 'n' lines in the input file
//
// Return 0 on success or -EIO on I/O error
//
int skiplines(FILE *fp, int n)
{
    int rc;
    int ch;

    while (n--) {
        if (feof(fp)) {
            break;
        }

        while ((ch = fgetc(fp)) != EOF && ch != '\r' && ch != '\n') {
        }

        if ((rc = skipeol(fp, ch)) != 0) {
            return rc;
        }
    }

    return 0;
}

// print the rest of the line the file pointer currently points at
//
int printline(FILE *fp)
{
    int rc;
    int ch;

    if (feof(fp)) {
        return 0;
    }

    while ((ch = fgetc(fp)) != EOF && ch != '\r' && ch != '\n') {
        fputc(ch, stdout);
    }
    fputc('\n', stdout);

    if ((rc = skipeol(fp, ch)) != 0) {
        return rc;
    }

    return 0;
}

// Compare the current lines starting at  fp[0] and fp[1].
// Return 0 on match and 1 on mismatch. Even if the lines
// don't match, both file pointers will be advanced past the
// terminating newline.
//
// Returns -EIO on I/O error on either file
//
int complines(FILE *fp[])
{
    int rc = 0;
    int skiprc;
    int ch[2];
    int i;

    for (;;) {
        ch[0] = fgetc(fp[0]);
        ch[1] = fgetc(fp[1]);
        
        if ((ch[0] == EOF && ferror(fp[0])) || (ch[1] == EOF && ferror(fp[1]))) {
            return -EIO;
        }

        if (ch[0] != ch[1]) {
            rc = 1;
        }

        if (ch[0] == EOF || ch[0] == '\r' || ch[0] == '\n') {
            break;
        }
    }

    for (i = 0; i < 2; i++) {
        while (ch[i] != EOF && ch[i] != '\n' && ch[i] != '\r') {
            ch[i] = fgetc(fp[i]);
        }

        if ((skiprc = skipeol(fp[i], ch[i])) != 0) {
            return skiprc;
        }
    }

    return rc;
}

//
// Parse the right-hand file into (serial, hash) pairs where serial is 
// the one-based line number.
//
// Returns 0 on success 
// -ENOMEM on out of memory
// -EIO on I/O error 
//
int parse_right_file(FILE *fp, vec_t *V)
{
    int rc;
    unsigned hash;
    serhash_t *sh;

    if ((rc = SH_INIT(V)) != 0) {
        return rc;
    }

    if ((sh = SH_APPEND(V)) == NULL) {
        return -ENOMEM;
    }
    sh->hash = ~1;
    sh->serial = 0;

    while ((rc = next_line(fp, &hash)) == 1) {
        if ((sh = SH_APPEND(V)) == NULL) {
            return -ENOMEM;
        }
        sh->serial = V->n - 1;
        sh->hash = hash;
    }

    if (rc < 0) {
        return rc;
    }

    return vec_clamp(V);
}

//
// Build the equivalence classes for a vector of (serial, hash) pairs.
//
// Return 0 on success or -ENOMEM on out of memory
//
int buildeq(vec_t *V, vec_t *E)
{
    int rc;
    int i;
    serhash_t *sh;
    eclass_t *ec;

    if ((rc = EV_INIT(E)) != 0) {
        return rc;
    }

    if ((ec = EV_APPEND(E)) == NULL) {
        return -ENOMEM;
    }

    ec->serial = 0;
    ec->last = 1;
    
    sh_sort(V);
    ec = EV_AT(E, 0);
    ec->last = 1;

    for (i = 1; i < V->n; i++) {
        sh = SH_AT(V, i);
        if ((ec = EV_APPEND(E)) == NULL) {
            return -ENOMEM;
        }
        ec->serial = sh->serial;
        ec->last = (i == V->n - 1) || (sh->hash != sh[1].hash); 
    }

    return vec_clamp(E);
}

//
// Parse the left hand file in terms of the right hand file's
// eq classes.
//
// Returns 0 on success 
// -ENOMEM on out of memory
// -EIO on I/O error 
//
int parse_left_file(FILE *fp, vec_t *V, vec_t *E, vec_t *P)
{
    int rc;
    unsigned hash;
    int *pi;
    int i;

    if ((rc = INTV_INIT(P)) != 0) {
        return rc;
    }

    if ((pi = INTV_APPEND(P)) == NULL) {
        return -ENOMEM;
    }
    *pi = 0;

    while ((rc = next_line(fp, &hash)) == 1) {
        i = sh_search(SH_AT(V, 0), EV_AT(E, 0), V->n, hash);
        
        if ((pi = INTV_APPEND(P)) == NULL) {
            return -ENOMEM;
        } 

        *pi = i;
    }

    return vec_clamp(P);
}

// Search between K[low] and K[high] for an index where 
// K[s]->b < j and K[s+1]->b > j. Note that this implies
// that K[high+1] must exist.
//
// Return the index s if found, else -1.
//
int merge_search(candidate_t **K, int low, int high, int j)
{
    // K is kept in order of increasing b so we can
    // use binary search 
    //
    while (low <= high)
    {
        int mid = (low + high) / 2;
        int rs = K[mid]->b;
        int re = K[mid+1]->b;
        if (j > rs && j < re) {
            return mid;
        } 

        if (j < rs) {
            high = mid - 1;
        } else {
            low = mid + 1;
        }
    }

    return -1;
}

// Given a line number in the left-hand file, look at all identical
// lines in the right-handle file in increasing order of line nuber
// (via the equivalence relation data computed earlier) and merge them
// into k-candidate chains. 
//
// K is array of rightmost K-candidates
// k is the highest k-candidate found
// i is the index into the left string
// E is the list of equivalence classes
// p is the index of the first eq class matching left[i] in E
//
// Returns 0 on success or -ENOMEM on out of memory
//
int merge(candidate_t **K, int *pk, int i, eclass_t *E, int p)
{
    int l = 0;

    candidate_t *c = K[0];
    candidate_t *t;
    int lowk = 0;
    int k = *pk;

    for (; !l; l = E[p].last, p++) {
        int j = E[p].serial;

        // criteria for a k-candidate
        // (1) Ai == Bj (we already know this is true)
        // (2) LCS exists between the A[:i] and B[:j]
        // (3) No common sequence of length k if either i or j reduced
        int fndk = merge_search(K, lowk, k, j);
        if (fndk == -1) {
            break;
        }

        // now, (i,j) is a fndk+1-candidate
        t = K[fndk];
        K[lowk] = c;        // save previous new entry
        c = candidate(i, j, t);
        if (c == NULL) {
            return -ENOMEM;
        }
        lowk = fndk + 1;

        if (fndk == k) {
            k++;
            K[k+1] = K[k];
            break;
        }
    }

    K[lowk] = c;
    *pk = k;
    return 0;
}

// Take the equivalence relations between the left and right files
// and compute the actual diff in the form of k-candidate chains
//
// Returns 0 on success or -ENOMEM on out of memory
//
int dodiff(vec_t *E, vec_t *P, candidate_t ***K, int *pk)
{
    int rc;
    int n = E->n;
    int m = P->n;
    int k_size = MIN(m, n) + 2;
    int k = 0;
    int i;

    if ((*K = (candidate_t **)malloc(k_size * sizeof(candidate_t *))) == NULL) {
        return -ENOMEM;
    }

    if (((*K)[0] = candidate(0, 0, NULL)) == NULL) {
        return -ENOMEM;
    }

    if (((*K)[1] = candidate(m + 1, n + 1, NULL)) == NULL) {
        return -ENOMEM;
    }

    int *pp = INTV_AT(P, 0);
    for (i = 1; i <= m; i++) {
        if (pp[i]) {
            rc = merge(*K, &k, i, EV_AT(E, 0), pp[i]);
            if (rc != 0) {
                return rc;
            }
        }
    }

    *pk = k;
    return 0;
}

// Scan for jackpots (dissimilar lines which have the same hash)
// If these are in J, then the lines aren't really equivalent and
// we remove the relation.
//
// Return 0 on success or -EIO on I/O error
//
int jackpot(int m, int *J, FILE *fp[])
{
    int rc;
    rewind(fp[0]);
    rewind(fp[1]);
    int l = 1;
    int r = 1;
    int i;

    for (i = 1; i <= m; i++) {
        if (J[i] == 0) {
            continue;
        }

        if (l < i) {
            if ((rc = skiplines(fp[0], i - l)) != 0) {
                return rc;
            }
            l = i;
        }

        if (r < J[i]) {
            if ((rc = skiplines(fp[1], J[i] - r)) != 0) {
                return rc;
            }
            r = J[i];
        }

        if ((rc = complines(fp)) < 0) {
            return rc;            
        }
        if (rc != 0) {
            printf("JACKPOT!\n");
            J[i] = 0;
        }

        l++, r++;
    }

    return 0;
}

// Take the longest candidate chain and transform it into an array,
// indexed by lines numbers of the left file, containing the matching
// line number in the right file, or 0 if there is no matching line. 
//
// Returns the match array or NULL on out of memory
//
int *build_matches(int m, int n, candidate_t *kk)
{
    int i;
    int *J;
    candidate_t *pc;

    J = (int*)malloc(sizeof(int*) * (m+1));
    if (J == NULL) {
        return J;
    }

    for (i = 0; i <= m; i++) {
        J[i] = 0;
    }

    for (pc = kk; pc != NULL; pc = pc->prev) {
        J[pc->a] = pc->b;
    }
    J[m] = n;

    return J;
}

// From the i <-> j array matching lines between the two files,
// build a representation of the differences in terms of blocks
// of lines in each file.
//
// Return 0 on success or -ENOMEM on out of memory
//
int build_deltas(vec_t *dv, int m, int *J)
{
    int rc;
    
    if ((rc = DV_INIT(dv)) != 0) {
        return rc;
    }

    int l = 0, r = 0;
    for (int i = 1; i <= m; i++) {
        delta_t *delta;

        if (J[i] == 0) {
            continue;
        }
        
        if (l+1 == i && r+1 == J[i]) {
            l++;
            r++;
            continue;
        }

        if ((delta = DV_APPEND(dv)) == NULL) {
            return -ENOMEM;
        }
        delta->left.start = l + 1;
        delta->left.end = i;
        delta->right.start = r + 1;
        delta->right.end = J[i];
        
        l = i;
        r = J[i];
    }

    return vec_clamp(dv);
}

// From the given file pointer, skip the next 'skip' lines, and
// the print the next 'print' lines to stdout, each prefaced with
// 'prefix' and a space.
//
// Return 0 on success or -EIO on I/O error
//
int printrange(FILE *fp, int skip, int print, char prefix)
{
    int rc;

    if ((rc = skiplines(fp, skip)) != 0) {
        return rc;
    }

    while (print--) {
        printf("%c ", prefix);
        if ((rc = printline(fp)) != 0) {
            return rc;
        }
    }

    return 0;
}

// Print the diff in ed-command format
//
// Return 0 on success or -EIO on I/O error
//
int printdiff(FILE *fp[], vec_t *dv)
{
    int rc;
    int l = 1;
    int r = 1;

    rewind(fp[0]);
    rewind(fp[1]);

    for (int i = 0; i < dv->n; i++) {
        delta_t *delta = DV_AT(dv, i);

        if (delta->right.start == delta->right.end) {
            printf("%d", delta->left.start);
            if (delta->left.end > delta->left.start + 1) {
                printf(",%d", delta->left.end - 1);
            }
            printf("d%d\n", delta->right.start - 1);
            rc = printrange(fp[0], delta->left.start - l, delta->left.end - delta->left.start, '<');
            if (rc != 0) {
                return rc;
            }
            l = delta->left.end;
        } else if (delta->left.start == delta->left.end) {
            printf("%da%d", delta->left.start - 1, delta->right.start);
            if (delta->right.end > delta->right.start + 1) {
                printf(",%d", delta->right.end - 1);
            }
            printf("\n"); 
            rc = printrange(fp[1], delta->right.start - r, delta->right.end - delta->right.start, '>');
            if (rc != 0) {
                return rc;
            }
            r = delta->right.end;
        } else {
            printf("%d", delta->left.start);
            if (delta->left.end > delta->left.start + 1) {
                printf(",%d", delta->left.end - 1);
            }
            printf("c%d", delta->right.start);
            if (delta->right.end > delta->right.start + 1) {
                printf(",%d", delta->right.end - 1);
            }
            printf("\n"); 
            rc = printrange(fp[0], delta->left.start - l, delta->left.end - delta->left.start, '<');
            if (rc != 0) {
                return rc;
            }
            l = delta->left.end;
            rc = printrange(fp[1], delta->right.start - r, delta->right.end - delta->right.start, '>');
            if (rc != 0) {
                return rc;
            }
            r = delta->right.end;
        }
    }

    return 0;
}
 
// Print the (minimal) usage format
//
void usage(void)
{
    fprintf(stderr, "diff: file-1 file-1\n");
}

// free all memory used by diff
//
void freediff(diffstate_t *ds)
{
    int i;
    for (i = 0; i < 2; i++) {
        if (ds->fp[i]) {
            fclose(ds->fp[i]);
        }
    }

    vec_free(&ds->V);
    vec_free(&ds->E);
    vec_free(&ds->P);
    free(ds->K);
    free(ds->J);
    vec_free(&ds->dv);
}

int main(int argc, char **argv)
{
    diffstate_t ds;
    int i;
    int m, n;
    int k;
    int rc;

    if (argc != 3) { 
        usage();
        return 1;
    }

    memset(&ds, 0, sizeof(ds));

    for (i = 0; i < 2; i++) {
        const char *fn = argv[i+1];
        ds.fp[i] = fopen(fn, "rb");
        if (!ds.fp[i]) {
            perror(fn);
            return 1;
        }
    } 

    if ((rc = parse_right_file(ds.fp[1], &ds.V)) != 0) {
        goto fail;
    }

    sh_sort(&ds.V);
    
    if ((rc = buildeq(&ds.V, &ds.E)) != 0) {
        goto fail;
    }
    
    if ((rc = parse_left_file(ds.fp[0], &ds.V, &ds.E, &ds.P)) != 0) {
        goto fail;
    }

    vec_free(&ds.V);

    n = ds.E.n; 
    m = ds.P.n;

    if ((rc = dodiff(&ds.E, &ds.P, &ds.K, &k)) != 0) {
        goto fail;
    }

    if ((ds.J = build_matches(m, n, ds.K[k])) == NULL) {
        rc = -ENOMEM;
        goto fail;
    }

    if ((rc = jackpot(m, ds.J, ds.fp)) != 0) {
        goto fail;
    }

    vec_t dv;
    if ((rc = build_deltas(&ds.dv, m, ds.J)) != 0) {
        goto fail;
    }

    if ((rc = printdiff(ds.fp, &ds.dv)) != 0) {
        goto fail;
    }

    rc = ds.dv.n == 0 ? 0 : 1;
    freediff(&ds);
    return rc;

fail:
    errno = -rc;
    perror("diff");
    return 2;
}
