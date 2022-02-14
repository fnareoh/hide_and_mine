#include "functions.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

unsigned int sut( char * x, int * C, int k, char hash, char * y )
{
    int n = strlen(x);
    int i, j, l;
    
    y[0]='\0'; j = n;
    
    for ( i = 0; i < n; i++ )
    {
        if ( C[i] == 0 )
        {
            j = i;
            break;
        }
    }
    l = 0;
    if ( j + k - 1 < n )
        while ( l < k ) y[l++] = x[j++];
    
    int f, p, c;
    while ( j < n )
    {
        p = j - k; c = p + 1;
        if ( C[p] == 0 && C[c] == 0 )
            y[l++]=x[j++];
        if ( C[p] == 1 && C[c] == 0 )
        {
            if ( strncmp ( &x[c], &x[f], k-1 ) == 0)	y[l++] = x[j++];
            else
            {
                y[l++] = hash;
                for ( i = c; i < c + k - 1; i++ )
                    y[l++] = x[i];
                y[l++] = x[j++];
            }
        }
        if ( C[p] == 1 && C[c] == 1 )	j++;
        if ( C[p] == 0 && C[c] == 1 )
        {
            j++;
            f = c;
        }
    }
    y[l] = '\0';
    
    return l;
}
