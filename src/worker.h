#pragma once
#include "util/SimParams.h"
#include "util/RRuntime.h"
#include "model/Grid.h"

using namespace std;
using namespace poet;
/*Functions*/
void worker_function(t_simparams *params);


/*Globals*/
#define TAG_WORK 42
#define TAG_FINISH 43
#define TAG_TIMING 44
#define TAG_DHT_PERF 45
#define TAG_DHT_STATS 46
#define TAG_DHT_STORE 47
