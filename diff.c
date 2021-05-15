#include <stdio.h>
#include <stdlib.h>

#ifndef MIN
# define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

typedef struct serhash_t serhash_t;
struct serhash_t {
    unsigned serial;
    unsigned hash;
};

typedef struct eclass_t eclass_t;
struct eclass_t {
    unsigned serial;
    int last;
};

typedef struct vec_t vec_t;
struct vec_t {
    size_t elesize;
    unsigned alloc;
    unsigned n;
    char *data;
};

typedef struct candidate_t candidate_t;
struct candidate_t {
    int a;
    int b;
    candidate_t *prev;
};

FILE *fp[2];

/*****************************************************************************
 * Memory routines
 */
/**
 * Immediately ditch if we run out of memory.
 */
void oom_abort()
{
    fprintf(stderr, "\nout of memory\n");
    exit(2);
}

/**
 * Allocate memory, guaranteed to succeed or exit.
 */
void *safe_malloc(size_t n)
{
    void *p = malloc(n);
    if (!p) {
        oom_abort();
    }
    return p;
}

/**
 * Reallocate memory, guaranteed to succeed or exit.
 */
void *safe_realloc(void *p, size_t n)
{
    p = realloc(p, n);
    if (!p) {
        oom_abort();
    }
    return p;
}

/**
 * Free an allocated pointer
 */
void safe_free(void *p)
{
    free(p);
}

/*****************************************************************************
 * Generic vector  
 */

/**
 * Initialize an empty vector
 */
void vec_init(vec_t *pv, size_t elesize)
{
    const unsigned INITIAL_ALLOC = 16;

    pv->elesize = elesize;
    pv->alloc = INITIAL_ALLOC;
    pv->n = 0;
    pv->data = safe_malloc(elesize * INITIAL_ALLOC);
}

/**
 * Add an element to a node vector
 */
void *vec_append(vec_t *pv)
{
    size_t oldsize, newsize;

    if (pv->n == pv->alloc) {
        oldsize = pv->alloc * pv->elesize;
        pv->alloc *= 2;
        newsize = pv->alloc * pv->elesize;
        if (newsize <= oldsize) {
            // file way too many lines
            oom_abort();
        }
        pv->data = safe_realloc(pv->data, newsize);
    }

    pv->n++;
    return &pv->data[pv->elesize * (pv->n - 1)];
}

/**
 * Free a vector's data
 */
void vec_free(vec_t *v)
{
    safe_free(v->data);
    v->data = NULL;
    v->alloc = 0;
    v->n = 0;
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
 * serial-hash functions
 */

/**
 * Compare two serial-hash nodes on hash, then serial
 */
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

/**
 * Sort a serial-hash vector on hash, then serial
 */
void sh_sort(vec_t *sh)
{
    // NB we skip the first placeholder node
    qsort(sh->data + sizeof(serhash_t), sh->n - 1, sizeof(serhash_t), &sh_compare);
}

/**
 * Search for the last occurance of an equivalence class (based on hash)
 * in the given sh vector. Return the index or 0 if there is no match.
 */
unsigned sh_search(serhash_t *V, eclass_t *E, int n, unsigned hash)
{
    unsigned low = 1;
    unsigned high = n - 1;

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

/**
 * Read the next line from 'fp' and compute its hash.
 * 
 * Return 1 on success with 'phash' containing the hash.
 * Return 0 on EOF w/ no data to hash. 
 */
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

/**
 * Parse the right-hand file into (serial, hash) pairs where serial is 
 * the one-based line number. 
 */
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
}

/**
 * Build the equivalence classes for a vector of (serial, hash) pairs.
 */
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
}

/**
 * Parse the left hand file in terms of the right hand file's
 * eq classes.
 */
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
}

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
        int fndk = -1;
        for (int ik = lowk; ik <= k; ik++) {
            if (K[ik]->b < j && K[ik+1]->b > j) {
                fndk = ik;
                break;
            }
        }

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

void usage(void)
{
    fprintf(stderr, "diff: file-1 file-1\n");
}

int main(int argc, char **argv)
{
    int i;
    vec_t V;
    vec_t E;
    vec_t P;
    int m, n;
    int k_size, k;
    candidate_t **K;
    int *pp;

    if (argc != 3) { 
        usage();
        return 1;
    }

    for (i = 0; i < 2; i++) {
        const char *fn = argv[i+1];
        fp[i] = fopen(fn, "rb");
        if (!fp[i]) {
            perror(fn);
            return 1;
        }
    } 

    parse_right_file(fp[1], &V);
    sh_sort(&V);
    buildeq(&V, &E);
    parse_left_file(fp[0], &V, &E, &P);

    vec_free(&V);

    n = E.n;
    m = P.n;

    k_size = MIN(m, n) + 2;
    K = (candidate_t **)safe_malloc(k_size * sizeof(candidate_t *));
    k = 0;
    K[0] = candidate(0, 0, NULL);
    K[1] = candidate(m + 1, n + 1, NULL);

    pp = INTV_AT(&P, 0);
    for (i = 1; i <= m; i++) {
        if (pp[i]) {
            k = merge(K, k, i, EV_AT(&E, 0), pp[i]);
        }
    }

    candidate_t *pr, *cu, *nx;

    pr = candidate(m, n, NULL);
    cu = K[k];

    while (cu) {
        nx = cu->prev;
        cu->prev = pr;
        pr = cu;
        cu = nx;
    }

    int l = -1, r = -1;

    for (; pr; pr = pr->prev) {
        if (l+1 != pr->a || r+1 != pr->b) {
            if (r+1 == pr->b) {
                printf("%d", l+1);
                if (l+2 < pr->a) {
                    printf(",%d", pr->a-1);
                }
                printf("d%d\n", r);
            } else if (l+1 == pr->a) {
                printf("%da%d", l, r+1);
                if (r+2 < pr->b) {
                    printf("%d", pr->b-1);
                }
                printf("\n");
            } else {                
                printf("%d", l+1);
                if (l+2 < pr->a) {
                    printf(",%d", pr->a-1);
                }
                printf("c");
                printf("%d", r+1);
                if (r+2 < pr->b) {
                    printf("%d", pr->b-1);
                }
                printf("\n");
            }
        }

        l = pr->a;
        r = pr->b;
    }
}
