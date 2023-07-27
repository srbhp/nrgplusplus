How to build a Model
================

``examples`` directory contains examples of different Models. 

- Copy a folder here from ``examples`` to start woking on a specific model.
-  Edit ``CMakeLists.txt`` to add the model folder.
    for example ``add_subdirectory(ExampleModel)``  to the end of the ``CMakeLists.txt`` file.
- `build` the project by invoking `make`.
   - create a build directory and invoke cmake from there
     - ``mkdir build ; cd build ; cmake .. ; make``
     - if we manage to compile the project, then the executable will be in the ``build/ExampleModel/`` directory.
