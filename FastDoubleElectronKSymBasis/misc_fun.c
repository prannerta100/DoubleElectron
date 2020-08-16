#include "helper_functions.h"

int parity(int x){
    if(x % 2 == 0) return 1;
    return -1;
}

int intmax(int a, int b){
    return (a>b)?a:b;
}

int intmin(int a, int b){
    return (a<b)?a:b;
}

int check(int *lst, int sz, int x){
    int i;
    for(i=0;i<sz;i++){
        if(lst[i]==x) return i;
    }
    return -1;
}
