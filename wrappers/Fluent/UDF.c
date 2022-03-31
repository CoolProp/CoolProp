#include "udf.h"

double Props(char*, char, double, char, double, char*);

DEFINE_ON_DEMAND(call_coolprop) {
    real p, t, r;

    p = 100000.0;
    t = 300.0;

    r = Props((char*)"D", 'T', t, 'P', p / 1000, (char*)"Air");

    Message("p = %lf, t = %lf => r = %lf\n", p, t, r);
}
