# GPOSolver
 Repository for GPOSolver, add support for user-define dynamic problem definination besides static one which official provided.

## Official Documents
http://cmp.felk.cvut.cz/gposolver/data/gposolver-1.2.pdf

## Recommend Dependency in Ubuntu
Using SDPA Solver:
```bash
sudo apt install libsdpa-dev
```

## Usage
If you want to define dynamic problem definination, here is the steps.

1. (Options) Generate C++ Files by Matlab scripts in `Matlab` folder.

2. Modify Function `void init(const int n)` in `include/gposolver/gpoproblem.h`
    - const int n: dynamic parameters defined by your own.
    - function content: use your dynamic parameters to generate other parameters.

3. Modify Function `void init_function(const int n)` in `include/gposolver/gposolver_base.h`
    - const int n: dynamic parameters defined by your own.
    - function content: `problem.init(n)` -> `problem.init({your parameters})`

4. Set solver in `solver.cpp` with your dynamic parameters initialization.

## Credits
[GpoSolver](http://cmp.felk.cvut.cz/gposolver/)
