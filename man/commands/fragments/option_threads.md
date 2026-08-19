`--threads` *positive integer*
: Set the number of computation threads to use, from 1 to 1024. The
  number of threads should not exceed the number of available CPU
  cores. The value 0 is also accepted and, like the default, uses all
  available cores; decimal values are truncated to their integer
  part. On Linux, "available" accounts for the CPU affinity mask and
  the cgroup CPU quota of the running process, so a job confined by
  **taskset**(1), Slurm, Docker or Kubernetes launches one thread per
  core it was actually granted, rather than one per core the machine
  has.
