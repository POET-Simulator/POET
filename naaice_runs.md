# NAAICE Measurements

### Barite

1. With AI Surrogate and single retraining
```NAA_SPEC=10.3.10.42:12345:1:5 mpirun -n 24 ./poet --ai --ai-backend=2 --fn=1 ./barite_200_rt.R ./barite_200.qs2 output_barite_ai```

2. With copy/expert knowledge functionality
```NAA_SPEC=10.3.10.42:12345:1:5 mpirun -n 24 ./poet -c ./barite_200_rt.R ./barite_200.qs2 output_barite_expert_knowledge```

### Dolomite

1. With AI Surrogate and single retraining
```NAA_SPEC=10.3.10.42:12345:2:5 mpirun -n 24 ./poet --ai --ai-backend=2 --fn=2 ./dolo_interp_rt_dt2000.R ./dolo_interp.qs2 output_dolo_ai```

2. With copy/expert knowledge functionality
```NAA_SPEC=10.3.10.42:12345:2:5 mpirun -n 24 ./poet -c ./dolo_interp_rt_dt2000.R ./dolo_interp.qs2 output_dolo_expert_knowledge```