#pragma once
#include <RInside.h>

using namespace std;
using namespace Rcpp;

/*Functions*/
void worker_function(RInside &R);


/*Globals*/
#define TAG_WORK 42
#define TAG_FINISH 43
#define TAG_TIMING 44
#define TAG_DHT_PERF 45
#define TAG_DHT_STATS 46
#define TAG_DHT_STORE 47