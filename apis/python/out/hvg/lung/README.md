# Census `lung` cell stats

## Memory stats

### 50k rows per Dask task, 16 processes: [tdbs_50k_d16x1_T4.json](tdbs_50k_d16x1_T4.json)
<!-- `bmdf mem.sh tdbs_50k_d16x1_T4.json` -->
```bash
mem.sh tdbs_50k_d16x1_T4.json
# scheduler 39G
# worker_0 7.0G
# worker_1 7.2G
# worker_2 7.3G
# worker_3 7.2G
# worker_4 7.2G
# worker_5 7.2G
# worker_6 8.3G
# worker_7 5.0G
# worker_8 7.3G
# worker_9 7.5G
# worker_10 11G
# worker_11 7.5G
# worker_12 7.2G
# worker_13 7.0G
# worker_14 7.0G
# worker_15 7.8G
```

### 100k rows per Dask task, 16 processes: [tdbs_100k_d16x1_T4.json](tdbs_100k_d16x1_T4.json)
<!-- `bmdf mem.sh tdbs_100k_d16x1_T4.json` -->
```bash
mem.sh tdbs_100k_d16x1_T4.json
# scheduler 39G
# worker_0 11G
# worker_1 9.8G
# worker_2 9.3G
# worker_3 11G
# worker_4 9.9G
# worker_5 11G
# worker_6 8.8G
# worker_7 11G
# worker_8 9.8G
# worker_9 8.8G
# worker_10 8.8G
# worker_11 9.8G
# worker_12 10G
# worker_13 9.7G
# worker_14 11G
# worker_15 9.3G
```


### 200k rows per Dask task, 16 processes: [tdbs_200k_d16x1_T4.json](tdbs_200k_d16x1_T4.json)
<!-- `bmdf mem.sh tdbs_200k_d16x1_T4.json` -->
```bash
mem.sh tdbs_200k_d16x1_T4.json
# scheduler 39G
# worker_0 11G
# worker_1 15G
# worker_2 10G
# worker_3 12G
# worker_4 9.0G
# worker_5 9.0G
# worker_6 9.8G
# worker_7 12G
# worker_8 11G
# worker_9 11G
# worker_10 9.0G
# worker_11 11G
# worker_12 11G
# worker_13 11G
# worker_14 11G
# worker_15 11G
```

