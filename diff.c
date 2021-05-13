#include <stdio.h>
#include <stdlib.h>

typedef struct eclass_t eclass_t;

struct eclass_t {
    unsigned serial;
    unsigned hash;
    int last;
};

typedef struct vec_t vec_t;

struct vec_t {
    size_t elesize;
    unsigned alloc;
    unsigned n;
    char *data;
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
    exit(1);
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
 * equivalence classes
 */

#define EV_INIT(ev) vec_init(ev, sizeof(eclass_t))
#define EV_APPEND(ev) ((eclass_t*)vec_append(ev))
#define EV_AT(ev, idx) \
    ( ((eclass_t*)((ev)->data)) + (idx) )

/**
 * Compare two eq class nodes on hash, then serial
 */
int ev_compare(const void *l, const void *r)
{
    const eclass_t *ln = (const eclass_t *)l;
    const eclass_t *rn = (const eclass_t *)r;

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

/*****************************************************************************
 * vector of integers
 */
#define INTV_INIT(iv) vec_init(iv, sizeof(int))
#define INTV_APPEND(iv) ((int*)(vec_append(iv)))
#define INTV_AT(iv, idx) \
    ( ((int*)((iv)->data)) + (idx) )

/**
 * Sort an eq class vector on hash, then serial
 */
void ev_sort(vec_t *ev)
{
    qsort(ev->data + sizeof(eclass_t), ev->n - 1, sizeof(eclass_t), &ev_compare);
}

/**
 * Search for the last occurance of an equivalence class (based on hash)
 * in the given ec vector. Return the index or 0 if there is no match.
 */
unsigned ev_search(vec_t *ev, unsigned hash)
{
    unsigned low = 1;
    unsigned high = ev->n - 1;

    while (low <= high) {
        unsigned mid = (low + high) / 2;
        unsigned mid_h = EV_AT(ev, mid)->hash;

        if (hash < mid_h) {
            high = mid - 1;
        } else if (hash > mid_h) {
            low = mid + 1;
        } else {
            while (!EV_AT(ev, mid-1)->last) {
                mid--;
            }

            return mid;
        }
    }

    return 0;
}

/**
 * Parse an input file into (serial, hash) tuples, where each tuple
 * represents a line and 'serial' is the 1-based line number.
 */
void parse_right_file(FILE *fp, vec_t *ev)
{
    const unsigned HASH0 = 5381;
    unsigned hash = HASH0;
    int ch;
    eclass_t *ec;

    EV_INIT(ev);
    ec = EV_APPEND(ev);
    ec->hash = ~1;
    ec->serial = 0;
    ec->last = 0;

    while ((ch = fgetc(fp)) != EOF) {
        if (ch != '\r' && ch != '\n') {
            hash = ((hash << 5) + hash) + ch;
            continue;
        }

        ec = EV_APPEND(ev);
        ec->serial = ev->n - 1;
        ec->hash = hash;
        ec->last = 0;

        hash = HASH0;

        if (ch == '\r') {
            ch = fgetc(fp);
            if (ch != '\n') {
                ungetc(ch, fp);
            }
        }
    }
}

/**
 * Build the equivalence classes for a vector of (serial, hash) pairs.
 * Sorts the vector in-place.
 */
void buildeq(vec_t *ev)
{
    int i;
    eclass_t *ec;

    ev_sort(ev);
    ec = EV_AT(ev, 0);
    ec->last = 1;

    for (i = 1; i < ev->n; i++) {
        ec = EV_AT(ev, i);
        ec->last = (i == ev->n - 1) || (ec->hash != ec[1].hash); 
    }
}

/**
 *
 */
void parse_left_file(FILE *fp, vec_t *ev, vec_t *iv)
{
    const unsigned HASH0 = 5381;
    unsigned hash = HASH0;
    int ch;
    int i, j;

    INTV_INIT(iv);
    *INTV_APPEND(iv) = 0;

    j = 1;
    while ((ch = fgetc(fp)) != EOF) {
        if (ch != '\r' && ch != '\n') {
            hash = ((hash << 5) + hash) + ch;
            continue;
        }

        i = ev_search(ev, hash);
        *INTV_APPEND(iv) = i;

        hash = HASH0;

        if (ch == '\r') {
            ch = fgetc(fp);
            if (ch != '\n') {
                ungetc(ch, fp);
            }
        }
    }
}

void usage(void)
{
    fprintf(stderr, "diff: file-1 file-1\n");
}

int main(int argc, char **argv)
{
    int i;
    vec_t V;
    vec_t P;

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
    buildeq(&V);
    parse_left_file(fp[0], &V, &P);
}