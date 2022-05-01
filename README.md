# TOP Project - Karman's whirlwind

---

### Last realized simulation

![Simulation animation](./simulation.gif)

--- 

### Quickstart

#### Build:

```bash
git clone https://github.com/byjtew/TOP-Project.git
cd TOP-Project
mkdir build
cd build
cmake .. -DCMAKE_C_COMPILER=`which mpicc`
make
```

#### Run:

`mpirun -np 4 ./bin/main ./bin/config`

**Generate simulation in a GIF** (if you have GNU-Parallel):

```bash
bash ./bin/parallel_gen_animate_gif_legacy.sh ./bin/resultat.raw ./output.gif
```

**Generate simulation in a GIF** (if you do not have GNU-Parallel):

```bash
sh ./bin/gen_animate_gif_legacy.sh ./bin/resultat.raw ./output.gif
```
