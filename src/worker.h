#pragma once
#include "util/RRuntime.h"

using namespace std;
using namespace poet;
/*Functions*/
void worker_function(RRuntime R);


/*Globals*/
#define TAG_WORK 42
#define TAG_FINISH 43
#define TAG_TIMING 44
#define TAG_DHT_PERF 45
#define TAG_DHT_STATS 46
#define TAG_DHT_STORE 47
