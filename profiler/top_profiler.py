import data

# Processes the set of previously written logs

threshold = 1.10  # Percent difference

db = data.FileBasedProfileDB()
actual_max_ts = 0
dt = db.find("python ann_data.py")
last_two = dt[-2:]
c = 0

for s in dt:
    new_db = sorted(dt, key=lambda ProfileData: ProfileData.timestamp)

    L = [1, 2]
    L[0] = dt[0].user_time_sec
    L[1] = dt[1].user_time_sec
    for i in range(0, len(dt)):
        print(f"{i} dt[{i}].user_time_sec = {dt[i].user_time_sec} ts {dt[i].timestamp}")
    print(f"L0 = {L[0]} L1 {L[1]}")

    if threshold * float(L[1]) < float(L[0]) or float(L[1]) > threshold * float(L[0]):
        raise SystemExit(f"Potential performance degradation detected {L[0]} va {L[1]}")
    print("No recent performance degradation detected")
