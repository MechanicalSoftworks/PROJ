/******************************************************************************
 * Project:  PROJ.4
 * Purpose:  This is primarily material originating from pj_obs_api.c
 *           (now proj_4D_api.c), that does not fit into the API
 *           category. Hence this pile of tubings and fittings for
 *           PROJ.4 internal plumbing.
 *
 * Author:   Thomas Knudsen,  thokn@sdfe.dk,  2017-07-05
 *
 ******************************************************************************
 * Copyright (c) 2016, 2017, 2018, Thomas Knudsen/SDFE
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *****************************************************************************/

#define FROM_PROJ_CPP

#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "geodesic.h"
#include "proj_internal.h"

#include "proj/internal/internal.hpp"

using namespace NS_PROJ::internal;

enum pj_io_units pj_left (PJ *P) {
    int u = P->inverted? P->right: P->left;
    if (u==PJ_IO_UNITS_CLASSIC)
        return PJ_IO_UNITS_PROJECTED;
    return static_cast<enum pj_io_units>(u);
}

enum pj_io_units pj_right (PJ *P) {
    int u = P->inverted? P->left: P->right;
    if (u==PJ_IO_UNITS_CLASSIC)
        return PJ_IO_UNITS_PROJECTED;
    return static_cast<enum pj_io_units>(u);
}



/**************************************************************************************/
int pj_has_inverse(PJ *P) {
/***************************************************************************************
Check if a a PJ has an inverse.
***************************************************************************************/
    return ( (P->inverted && (P->host->fwd || P->host->fwd3d || P->host->fwd4d) ) ||
             ( P->host->inv || P->host->inv3d || P->host->inv4d) );
}


/* Move P to a new context - or to the default context if 0 is specified */
void proj_context_set (PJ *P, PJ_CONTEXT *ctx) {
    if (nullptr==ctx)
        ctx = pj_get_default_ctx ();
    proj_assign_context (P, ctx);
}


void proj_context_inherit (PJ *parent, PJ *child) {
    if (nullptr==parent)
        proj_assign_context (child, pj_get_default_ctx());
    else
        proj_assign_context (child, pj_get_ctx(parent));
}



/*****************************************************************************/
char *pj_chomp (char *c) {
/******************************************************************************
Strip pre- and postfix whitespace. Inline comments (indicated by '#') are
considered whitespace.
******************************************************************************/
    size_t i, n;
    char *comment;
    char *start = c;

    if (nullptr==c)
        return nullptr;

    comment = strchr (c, '#');
    if (comment)
        *comment = 0;

    n = strlen (c);
    if (0==n)
        return c;

    /* Eliminate postfix whitespace */
    for (i = n - 1;  (i > 0) && (isspace (c[i]) || ';'==c[i]);  i--)
        c[i] = 0;

    /* Find start of non-whitespace */
    while (0 != *start  &&  (';'==*start  ||  isspace (*start)))
        start++;

    n = strlen (start);
    if (0==n) {
        c[0] = 0;
        return c;
    }

    memmove (c, start, n + 1);
    return c;
}



/*****************************************************************************/
char *pj_shrink (char *c) {
/******************************************************************************
Collapse repeated whitespace. Remove '+' and ';'. Make ',' and '=' greedy,
consuming their surrounding whitespace.
******************************************************************************/
    size_t i, j, n;

    /* Flag showing that a whitespace (ws) has been written after last non-ws */
    bool ws = false;

    if (nullptr==c)
       return nullptr;

    pj_chomp (c);
    n = strlen (c);
    if (n==0)
        return c;

    /* First collapse repeated whitespace (including +/;) */
    i = 0;
    bool in_string = false;
    for (j = 0;  j < n;  j++) {

        if( in_string ) {
            if( c[j] == '"' && c[j+1] == '"' ) {
                c[i++] = c[j];
                j++;
            } else if( c[j] == '"' ) {
                in_string = false;
            }
            c[i++] = c[j];
            continue;
        }

        /* Eliminate prefix '+', only if preceded by whitespace */
        /* (i.e. keep it in 1.23e+08) */
        if ((i > 0) && ('+'==c[j]) && ws)
            c[j] = ' ';
        if ((i==0) && ('+'==c[j]))
            c[j] = ' ';

        // Detect a string beginning after '='
        if( c[j] == '"' && i > 0 && c[i-1] == '=' ) {
            in_string = true;
            ws = false;
            c[i++] = c[j];
            continue;
        }

        if (isspace (c[j]) || ';'==c[j]) {
            if (false==ws && (i > 0))
                c[i++] = ' ';
            ws = true;
            continue;
        }
        else {
            ws = false;
            c[i++] = c[j];
        }
    }
    c[i] = 0;
    n = strlen(c);

    /* Then make ',' and '=' greedy */
    i = 0;
    for (j = 0;  j < n;  j++) {
        if (i==0) {
            c[i++] = c[j];
            continue;
        }

        /* Skip space before '='/',' */
        if ('='==c[j] || ','==c[j]) {
            if (c[i - 1]==' ')
               c[i - 1] = c[j];
            else
                c[i++] = c[j];
            continue;
        }

        if (' '==c[j] && ('='==c[i - 1] || ','==c[i - 1]) )
            continue;

        c[i++] = c[j];
    }
    c[i] = 0;
    return c;
}



/*****************************************************************************/
size_t pj_trim_argc (char *args) {
/******************************************************************************
Trim all unnecessary whitespace (and non-essential syntactic tokens) from the
argument string, args, and count its number of elements.
******************************************************************************/
    size_t i, m, n;
    pj_shrink (args);
    n = strlen (args);
    if (n==0)
        return 0;
    bool in_string = false;
    for (i = m = 0;  i < n;  i++) {
        if (in_string ) {
            if( args[i] == '"' && args[i+1] == '"' ) {
                i++;
            } else if( args[i] == '"' ) {
                in_string = false;
            }
        }
        else if (args[i] == '=' && args[i+1] == '"' ) {
            i++;
            in_string = true;
        }
        else if (' '==args[i]) {
            args[i] = 0;
            m++;
        }
    }
    return m + 1;
}



/*****************************************************************************/
char **pj_trim_argv (size_t argc, char *args) {
/******************************************************************************
Create an argv-style array from elements placed in the argument string, args.

args is a trimmed string as returned by pj_trim_argc(), and argc is the number
of trimmed strings found (i.e. the return value of pj_trim_args()). Hence,
    int argc    = pj_trim_argc (args);
    char **argv = pj_trim_argv (argc, args);
will produce a classic style (argc, argv) pair from a string of whitespace
separated args. No new memory is allocated for storing the individual args
(they stay in the args string), but for the pointers to the args a new array
is allocated and returned.

It is the duty of the caller to free this array.
******************************************************************************/

    if (nullptr==args)
        return nullptr;
    if (0==argc)
        return nullptr;


    /* turn the input string into an array of strings */
    char** argv = (char **) calloc (argc, sizeof (char *));
    if (nullptr==argv)
        return nullptr;
    for(size_t i = 0, j = 0; j < argc; j++) {
        argv[j] = args + i;
        char* str = argv[j];
        size_t nLen = strlen(str);
        i += nLen + 1;
    }
    return argv;
}


/*****************************************************************************/
std::string pj_double_quote_string_param_if_needed(const std::string& str) {
/*****************************************************************************/
    if( str.find(' ') == std::string::npos ) {
        return str;
    }
    return '"' + replaceAll(str, "\"", "\"\"") + '"';
}

/*****************************************************************************/
char *pj_make_args (size_t argc, char **argv) {
/******************************************************************************
pj_make_args is the inverse of the pj_trim_argc/pj_trim_argv combo: It
converts free format command line input to something proj_create can consume.

Allocates, and returns, an array of char, large enough to hold a whitespace
separated copy of the args in argv. It is the duty of the caller to free this
array.
******************************************************************************/

    try
    {
        std::string s;
        for( size_t i = 0; i < argc; i++ )
        {
            const char* equal = strchr(argv[i], '=');
            if( equal ) {
                s += std::string(argv[i], equal - argv[i] + 1);
                s += pj_double_quote_string_param_if_needed(equal + 1);
            } else {
                s += argv[i];
            }
            s += ' ';
        }

        char* p = pj_strdup(s.c_str());
        return pj_shrink (p);
    }
    catch( const std::exception& ) {
        return nullptr;
    }
}
