//
// TODO 
//   better memory handling for library
//   split into library interface
//
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
    unsigned alloc;         // # of elements allocated
    unsigned n;             // # of elements in use
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
 * Memory routines
 */

//
// Immediately ditch if we run out of memory.
//
void oom_abort()
{
    fprintf(stderr, "\nout of memory\n");
    exit(2);
}

//
// Allocate memory, guaranteed to succeed or exit.
//
void *safe_malloc(size_t n)
{
    void *p = malloc(n);
    if (!p) {
        oom_abort();
    }
    return p;
}

//
// Reallocate memory, guaranteed to succeed or exit.
//
void *safe_realloc(void *p, size_t n)
{
    p = realloc(p, n);
    if (!p) {
        oom_abort();
    }
    return p;
}

//
// Free an allocated pointer
//
void safe_free(void *p)
{
    free(p);
}

/*****************************************************************************
 * Generic vector  
 */

//
// Initialize an empty vector
//
void vec_init(vec_t *pv, size_t elesize)
{
    const unsigned INITIAL_ALLOC = 16;

    pv->elesize = elesize;
    pv->alloc = INITIAL_ALLOC;
    pv->n = 0;
    pv->data = safe_malloc(elesize * INITIAL_ALLOC);
}

//
// Add an element to a vector
//
void *vec_append(vec_t *pv)
{
    size_t oldsize, newsize;

    if (pv->n == pv->alloc) {
        oldsize = pv->alloc * pv->elesize;
        pv->alloc *= 2;
        newsize = pv->alloc * pv->elesize;
        if (newsize <= oldsize) {
            // size has wrapped around, too big!
            oom_abort();
        }
        pv->data = safe_realloc(pv->data, newsize);
    }

    pv->n++;
    return &pv->data[pv->elesize * (pv->n - 1)];
}

// clamp a vector so its allocated size is the exactly the number of elements
// in use. used to recover space reserved for new elements when it's known
// that the vector will never again need to grow.
//
void vec_clamp(vec_t *pv)
{
    pv->data = safe_realloc(pv->data, pv->elesize * pv->n);
    pv->alloc = pv->n;
}

//
// Free a vector's data
//
void vec_free(vec_t *v)
{
    if (v->data) {
        safe_free(v->data);
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
    candidate_t *p = (candidate_t *)safe_malloc(sizeof(candidate_t));
    p->a = a;
    p->b = b;
    p->prev = prev;
    return p;
}

/*****************************************************************************
 * file parsing and algorithm
 */

//
// Read the next line from 'fp' and compute its hash.
// 
// Return 1 on success with 'phash' containing the hash.
// Return 0 on EOF w/ no data to hash. 
//
int next_line(FILE *fp, unsigned *phash)
{
    int ch;
    const unsigned HASH0 = 5381;
    unsigned hash = HASH0;

    if (feof(fp)) {
        return 0;
    }

    while ((ch = fgetc(fp)) != EOF) {
        if (ch == '\r' || ch == '\n') {
            break;
        }
        hash = ((hash << 5) + hash) + ch;
    }

    // we read a blank link that ended in EOF - that counts as EOF
    //
    if (ch == EOF && hash == HASH0) {
        return 0;
    }

    if (ch == '\r') {
        ch = fgetc(fp);
        if (ch != '\n') {
            ungetc(ch, fp);
        }
    }

    *phash = hash;
    return 1;
}

// Skip the next 'n' lines in the input file
//
void skiplines(FILE *fp, int n)
{
    int ch;

    while (n--) {
        if (feof(fp)) {
            return;
        }

        while ((ch = fgetc(fp)) != EOF) {
            if (ch == '\r' || ch == '\n') {
                break;
            }
        }

        if (ch == '\r') {
            ch == fgetc(fp);
            if (ch != '\n') {
                ungetc(ch, fp);
            }
        }
    }
}

// print the rest of the line the file pointer currently points at
//
void printline(FILE *fp)
{
    int ch;

    if (feof(fp)) {
        return;
    }

    while ((ch = fgetc(fp)) != EOF) {
        if (ch == '\r' || ch == '\n') {
            break;
        }
        fputc(ch, stdout);
    }
    fputc('\n', stdout);

    if (ch == '\r') {
        ch == fgetc(fp);
        if (ch != '\n') {
            ungetc(ch, fp);
        }
    }
}

// Compare the current lines starting at  fp[0] and fp[1].
// Return 0 on match and non-zero on mismatch. Even if the lines
// don't match, both file pointers will be advanced past the
// terminating newline.
//
int complines(FILE *fp[])
{
    int rc = 0;
    int ch[2];
    int i;

    for (;;) {
        ch[0] = fgetc(fp[0]);
        ch[1] = fgetc(fp[1]);

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

        if (ch[i] == '\r') {
            ch[i] = fgetc(fp[i]);
            if (ch[i] != '\n') {
                ungetc(ch[i], fp[i]);
            }
        }
    }

    return rc;
}

//
// Parse the right-hand file into (serial, hash) pairs where serial is 
// the one-based line number. 
//
void parse_right_file(FILE *fp, vec_t *V)
{
    unsigned hash;
    serhash_t *sh;

    SH_INIT(V);
    sh = SH_APPEND(V);
    sh->hash = ~1;
    sh->serial = 0;

    while (next_line(fp, &hash)) {
        sh = SH_APPEND(V);
        sh->serial = V->n - 1;
        sh->hash = hash;
    }

    vec_clamp(V);
}

//
// Build the equivalence classes for a vector of (serial, hash) pairs.
//
void buildeq(vec_t *V, vec_t *E)
{
    int i;
    serhash_t *sh;
    eclass_t *ec;

    EV_INIT(E);
    ec = EV_APPEND(E);
    ec->serial = 0;
    ec->last = 1;
    
    sh_sort(V);
    ec = EV_AT(E, 0);
    ec->last = 1;

    for (i = 1; i < V->n; i++) {
        sh = SH_AT(V, i);
        ec = EV_APPEND(E);
        ec->serial = sh->serial;
        ec->last = (i == V->n - 1) || (sh->hash != sh[1].hash); 
    }

    vec_clamp(E);
}

//
// Parse the left hand file in terms of the right hand file's
// eq classes.
//
void parse_left_file(FILE *fp, vec_t *V, vec_t *E, vec_t *P)
{
    unsigned hash;
    int i;

    INTV_INIT(P);
    *INTV_APPEND(P) = 0;

    while (next_line(fp, &hash)) {
        i = sh_search(SH_AT(V, 0), EV_AT(E, 0), V->n, hash);
        *INTV_APPEND(P) = i;
    }

    vec_clamp(P);
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
int merge(candidate_t **K, int k, int i, eclass_t *E, int p)
{
    int l = 0;

    candidate_t *c = K[0];
    candidate_t *t;
    int lowk = 0;

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
        lowk = fndk + 1;

        if (fndk == k) {
            k++;
            K[k+1] = K[k];
            break;
        }
    }

    K[lowk] = c;

    return k;
}

// Take the equivalence relations between the left and right files
// and compute the actual diff in the form of k-candidate chains
//
void dodiff(vec_t *E, vec_t *P, candidate_t ***K, int *pk)
{
    int n = E->n;
    int m = P->n;
    int k_size = MIN(m, n) + 2;
    int k = 0;
    int i;

    *K = (candidate_t **)safe_malloc(k_size * sizeof(candidate_t *));

    (*K)[0] = candidate(0, 0, NULL);
    (*K)[1] = candidate(m + 1, n + 1, NULL);

    int *pp = INTV_AT(P, 0);
    for (i = 1; i < m; i++) {
        if (pp[i]) {
            k = merge(*K, k, i, EV_AT(E, 0), pp[i]);
        }
    }

    *pk = k;
}

// Scan for jackpots (dissimilar lines which have the same hash)
// If these are in J, then the lines aren't really equivalent and
// we remove the relation.
//
void jackpot(int m, int *J, FILE *fp[])
{
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
            skiplines(fp[0], i - l);
            l = i;
        }

        if (r < J[i]) {
            skiplines(fp[1], J[i] - r);
            r = J[i];
        }

        if (complines(fp) != 0) {
            printf("JACKPOT!\n");
            J[i] = 0;
        }

        l++, r++;
    }
}

// Take the longest candidate chain and transform it into an array,
// indexed by lines numbers of the left file, containing the matching
// line number in the right file, or 0 if there is no matching line. 
//
int *build_matches(int m, int n, candidate_t *kk)
{
    int i;
    int *J;
    candidate_t *pc;

    J = (int*)safe_malloc(sizeof(int*) * (m+1));
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
void build_deltas(vec_t *dv, int m, int *J)
{
    DV_INIT(dv);

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

        delta = DV_APPEND(dv);
        delta->left.start = l + 1;
        delta->left.end = i;
        delta->right.start = r + 1;
        delta->right.end = J[i];
        
        l = i;
        r = J[i];
    }

    vec_clamp(dv);
}

// From the given file pointer, skip the next 'skip' lines, and
// the print the next 'print' lines to stdout, each prefaced with
// 'prefix' and a space.
//
void printrange(FILE *fp, int skip, int print, char prefix)
{
    skiplines(fp, skip);
    while (print--) {
        printf("%c ", prefix);
        printline(fp);
    }
}

// Print the diff in ed-command format
//
void printdiff(FILE *fp[], vec_t *dv)
{
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
            printrange(fp[0], delta->left.start - l, delta->left.end - delta->left.start, '<');
            l = delta->left.end;
        } else if (delta->left.start == delta->left.end) {
            printf("%da%d", delta->left.start - 1, delta->right.start);
            if (delta->right.end > delta->right.start + 1) {
                printf(",%d", delta->right.end - 1);
            }
            printf("\n"); 
            printrange(fp[1], delta->right.start - r, delta->right.end - delta->right.start, '>');
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
            printrange(fp[0], delta->left.start - l, delta->left.end - delta->left.start, '<');
            l = delta->left.end;
            printrange(fp[1], delta->right.start - r, delta->right.end - delta->right.start, '>');
            r = delta->right.end;
        }
    }
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
    safe_free(ds->K);
    safe_free(ds->J);
    vec_free(&ds->dv);
}

int main(int argc, char **argv)
{
    diffstate_t ds;
    int i;
    int m, n;
    int k;

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

    parse_right_file(ds.fp[1], &ds.V);
    sh_sort(&ds.V);
    buildeq(&ds.V, &ds.E);
    parse_left_file(ds.fp[0], &ds.V, &ds.E, &ds.P);
    vec_free(&ds.V);

    n = ds.E.n; 
    m = ds.P.n;

    dodiff(&ds.E, &ds.P, &ds.K, &k);
    ds.J = build_matches(m, n, ds.K[k]);
    jackpot(m, ds.J, ds.fp);

    vec_t dv;
    build_deltas(&ds.dv, m, ds.J);

    printdiff(ds.fp, &ds.dv);

    freediff(&ds);
}
