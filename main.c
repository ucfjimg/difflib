#include "diff.h"

#include <errno.h>
#include <stdio.h>

// Print the (minimal) usage format
//
void usage(void)
{
    fprintf(stderr, "diff: file-1 file-1\n");
}



int main(int argc, char **argv)
{
    int i;
    diff_t d;

    if (argc != 3) { 
        usage();
        return 1;
    }

    for (i = 0; i < 2; i++) {
        const char *fn = argv[i+1];
        d.fp[i] = fopen(fn, "rb");
        if (!d.fp[i]) {
            perror(fn);
            return 1;
        }
    }

    i = diff(&d);
    if (i != 0) {
        errno = -i;
        perror("diff");
        return 2;
    }

    if (printdiff(&d) < 0) {
        return 2;
    }

    return d.ndelta == 0 ? 0 : 1;
}
