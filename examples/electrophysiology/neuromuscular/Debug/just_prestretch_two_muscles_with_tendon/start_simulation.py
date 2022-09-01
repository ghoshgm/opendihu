import subprocess
import shlex
import timeit

start = timeit.default_timer()



iterations = 1

for _ in range(iterations):


    # run muscle-simulation
    sim_endtime = 20

    sim_process_input = shlex.split(f"./prestretch_two_muscles_with_tendon ../settings.py --tend={sim_endtime}")

    subprocess.run(sim_process_input)

stop = timeit.default_timer()

print('Time: ', stop - start) 
