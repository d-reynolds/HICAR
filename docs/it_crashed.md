# What to do if the model crashes

Seemingly random or intractible model crashes can be one of the most frustrating parts of working with geophysical models. Here, we recommend some steps to take if you find yourself visiting this page.

## Inspect the logs

If you haven't done so already, inspect the model output, or any logs (.out or .err redirect files). If the model either gave a descriptive reason for the crash in .out, it may be a Feature Not A Bug™️, and the model is trying to stop you from a bad physics combination, or crashed due to a known reason. If not, .err file may also contain a clue already as to where in the code, or why, the crash occured. The [errors](errors.md) section includes some common reasons for a crash.

## Check your branch

If you are compiling code from the main branch, that code is intended to work, and you may be experiencing a real, unknown bug. If you would like to dig deeper and contribute to this open-source project, you can try to debug the error by following the steps herein. At the least, we would greatly appreciate it if you filed a bug report under issues so that the core developers are aware of it. The steps below will help you to generate more informative diagnostics to assist the developers in addressing your bug. 

If you are compiling code from a feature branch (or `main` between releases), hop to it soldier. Unreleased code may be new and shiny, but is not guaranteed to work. You may first try running some of the [tests](testing.md) shipped with the model, and a failure here may already be informative about what is breaking. In any case, this guide is intended to help you get to the bottom of things.

## Check your input data

If the model was running fine and crashed in the physics loop, check that your input data around that time is not corrupted in some way. If this is your first run with a new initial conditions file, check that the variables of that file are all reasonable. The model will try to perform some of these steps for you, but it ain't perfect.

## Make a small, reproducible test case

If your crash occured during initialization, it should be relatively easy to reproduce without a lengthy lead-up time for subsequent runs, so skip to the next section. If the crash occured after initialization though, take the following steps:

1. Output a restart file just before the model crash occurs
2. Set the model's `restart_date` to the date from the restart file, and set `restart = .True.` in the namelist
3. Rerun and check the outcome:

    a. The model crashes at the same time: **Congratulations, your crash is reproducible 🎉**

    **-- or --**

    b. The model crashed before/after the initial crash ➡️ your crash is not reproducible, continue to [non-reproducible crashes](#non-reproducible-crashes).

## Debug run

Assuming that you were using the default, release exe of HICAR, we will not go through the steps to build a debug exe. This will give us more information about the crash if it is caused by an array bounds error, or some other tractible error signal.

1. First create a debug build of HICAR

```bash
cd /path/to/hicar/repo
mkdir build_debug
cd build_debug
cmake ../ -DMODE=debug #Plus any other configuration flags you used in your crashing release build
make -j
make install
```

This will install a HICAR_debug (or HICAR_debug_gpu) executable to your bin/ directory

2. Set debug = .True. in the namelist. This will slow down the model signifficantly, but will trigger additional runtime checks after each step of the physics loop to catch any bad values right when they occur. If setting debug = .True. causes the run to progress far too slowly to be useful, you can still run with debug = .False.

3. Now rerun your small, reproducible test case. If you're lucky, the standard error output will now include more information about why the model crashed, and where in the code it crashed. If not, or if that location is deep inside of an archane physics library:

## Manual debug

This is where things become more art than science. At this stage, it is a good idea to file a github issue as well so that others are aware of the bug, and won't feel alone if they hit the same crash.

When we get into manual debugging, it is important to be rigorously honest with yourself about what you are testing. Our goal here is to have two runs, one run that produces an error, and another that does not, separated by only a single change in the model configuration, code, or namelist options. It can be tempting to work fast and change numerous variables at once, seeking to clear the bug/crash-site, but a passing run then tells us a little about what caused the crash. Conversely, iterating over small changes when you already have a hunch about what option or physics module is causing the crash will waste time. The best method is to start with larger steps to isolate a set of changes that clears the bug, and then pair down this set using additional testing, until a single change is identified which clears the bug. 

Atmospheric models tightly interleave various physical processes, so a bug in some part of the code may not trigger anything, appear physical on its first iteration, and later tank another part of the code. An array out-of-bounds error in a well-tested physics submodule, like a microphysics scheme, is a typical symptom of this kind of problem. Most of the physics modules in the model have been developed by independent groups and used in a variety of atmospheric models -- [Occam's Razor](https://en.wikipedia.org/wiki/Occam%27s_razor) says that the bug is most likely not in their code. The first steps we usually try are turning off various physics options. Since we will not be editing the namelist, save a copy of the original namelist that caused the crash. At the end of debugging, we will re-run with this namelist to prove to ourselves that the crash is gone. The most heinously intertwined parameterizations are the surface physics (Land Surface Model (LSM), Surface Layer (SFC), Water, Snow Model (SM), and PBL scheme). We will often start by turning off all of these, by setting their namelist options to "none" in the `&physics` namelist section:

```bash
    pbl = 'none'
    lsm = 'none'
    sfc = 'none'
    water = 'none'
    sm = 'none'
```

Try running again, and see if the model runs through your previous time of the crash. If it does, inspect the model output after the time where it crashed before. Do all of the fields look OK? Temperature, wind components? If there are NaN values here, or if the model still crashed, the issue is still there. The next stage of escalation is to turn off the radiation and microphysics options. You can do this one-by-one, or both together.

```bash
    rad = 'none'
    mp = 'none'
```

If you were running with `terrain_shading=.True.` in the `&rad` namelist section, you will now also need to set `terrain_shading=.False.`. Follow the same conditions outlined above: is the crash gone, and the model output looks clean? If not, the last set of options is to turn off the wind solver and advection, by setting:


```bash
    wind = 'none'
    adv = 'none'
```

It is also a good idea to set the options in the `&wind` section to their default values now as well. If a bug still persists, it is deeper in the model infrastructure code. These can be a pain to debug, and we recommend getting in touch with the core developers if you haven't already.

If one of the proceeding steps _did_ disable the bug though, you are closer to surrounding it. The next step would be changing any of the fine-grained namelist options related to the breaking physics module, for example disabling terrain shading if you found the bug could be toggled by setting `rad = 'none'`. 

These steps are not meant to be exhaustive, but are intended to guide new model users or programmers through some of the steps that we have found most helpful when debugging.

## Non-reproducible Crashes

If the crash is not reproducible using the steps outlined earlier, this suggests an error in I/O, communications, or your own hardware setup. These can either be ephemeral (in the case of limited hardware resources caused by multiple users on an unscheduled system), or extremely intractible (in the case of communications race conditions). We hope it's the first one too, but in any case file an issue and we can work on this together 🌈.