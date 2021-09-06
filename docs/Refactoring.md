# Patterns

## Restoring Data from File
An often-used trick to avoid having many function parameters is to store variables from the environment to disk and to restore them again later (e.g [here](https://github.com/RECETOX/recetox-xMSannotator/blob/c7ad3bb2f4e7cc6aa35d8a0804fe6ef7a8389729/R/multilevelannotationstep3.R#L14-L17)). This is also used to avoid having to export variables for parallel function calls.

### High Level Functions
Instead of reading the data from disk inside the high-level function, all variables stored in the file should be instead passed as function parameters, see [here](https://github.com/hechth/recetox-xMSannotator/blob/e9ea38a1654993488af8b188c2a39f008ecbbfd5/R/multilevelannotationstep2.R#L108-132).

### Parallel Function Calls
Parameters can be exported to the computing cluster as seen [here](https://github.com/hechth/recetox-xMSannotator/blob/e9ea38a1654993488af8b188c2a39f008ecbbfd5/R/multilevelannotation.R#L1264-L1274). It has to be made sure that the variables are exported from the correct environment, while the default is the parent environment.

## Parallelism
The code uses different parallel programming solutions, such as [`foreach`](https://www.rdocumentation.org/packages/foreach/versions/1.5.1) coupled with `%dopar%`, or `snow` with `parLapply(...)`. Both systems have advantages and disadvantages, but ideally, there should be a consensus on which system to use. Note that `snow` is outdated and superseded by `parallel` ([source](https://bookdown.org/rdpeng/rprogdatascience/parallel-computation.html)).

### `snow` & `parallel`
The `snow` and `parallel` packages can perform process communication via TCP, making them platform-independent. Variables can be passed to functions, but it seems that they have to be made available on the other processes, either via explicit `clusterExport(...)` function call, or by storing them in the global environment - this is my assumption, which hasn't been properly validated yet. Storing any variables in the global environment is **HIGHLY DISCOURAGED** as it messes with tests.

The `parallel` package provides a wrapper around multiple backends, `snow` being one of them.

### `foreach`
The `foreach` package provides a %dopar% special operator (see more on custom operators [here](https://stackoverflow.com/questions/25179457/r-what-are-operators-like-in-called-and-how-can-i-learn-about-them)) which allows a bit easier syntax. In order to use the parallel version, also a backend has to be registered (`snow` or `parallel`), otherwise it falls back to a sequential version, making the code a bit more robust.

The `foreach` operator provides options to combine the results, but the behaviour of the functions is not documented and produces somewhat bizarre results, see the [Possible Issues](https://github.com/RECETOX/recetox-xMSannotator/wiki/Possible-Issues) section on this wiki.

#### Debugging
Many of the `foreach` constructs call the code inline, meaning that there is no explicit function definition. This causes errors while setting breakpoints, therefore the code should be extracted into a function. Another issue with debugging is that the %dopar% operator causes a namespace related error when being used in the `Debug R-Package` configuration of `vscDebugger`. Before filing an issue to [vscDebugger](https://github.com/ManuelHentschel/VSCode-R-Debugger) or [foreach](https://github.com/RevolutionAnalytics/foreach) there should be consistent testing, ideally also whether this bug also appears in RStudio.