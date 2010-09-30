=====================================
Cookbook for Building TEAL Interfaces
=====================================

.. abstract::
   :author: Warren J. Hack, Chris Sontag
   :date: 23 September 2010
   
   The release of the Task Editor And Launcher(TEAL) with STScI_Python
   v2.10 in June 2010 provided the tools for building powerful GUI
   interfaces for editing the parameters of complex tasks and running that
   task with minimal effort. Learning how to use something new always
   takes a special effort, and this document provides a step-by-step
   walkthrough of how to build TEAL interfaces for any Python task to
   make this effort as easy as possible.
   
------------
Introduction
------------

The new TEAL GUI can be added to nearly any Python task that allows users to set parameters to control the operation of the task. Adding a TEAL interface to a Python task requires some minor updates to the task's code in order to allow TEAL to create and control the GUI for setting all the necessary parameters. TEAL itself relies on the `ConfigObj module`_ for the basic parameter handling functions, with additional commands for implementing enhanced logic for controlling the GUI itself based on parameter values. The GUI not only guides the user in setting the parameters, but also provides the capability to load and save parameter sets and the ability to read a help file while still editing the parameters.  The interface to TEAL can also be set up alongside a command-line interface to the task.  This document provides the basic information necessary for implementing a TEAL interface for nearly any Python task to take full advantage of the control it provides the user in setting the task parameters. 

This document does not assume the user has any familiarity with using configobj in any manner and as as result includes very basic information which developers with some experience with configobj can simply skip over. 

The development of the TEAL interface for the task `resetbits` in the `betadrizzle` package.  More elaborate examples will be explained after the development of the TEAL interface for `resetbits` has been described. 

----------------------
Building the Interface
----------------------

The order of operations provided by this document is not the only order in which these steps can be performed.  This order starts with the simplest operation then leads the developer into what needs to be done next with the least amount of iteration.  


Step 1: Defining the Parameters
===============================

The primary purpose for developing a TEAL interface is to provide a GUI which can be used to set the values for the task's parameters. This requires that the developer identify the full set of task parameters which the user will be required to provide when running the task. The signature for the task `resetbits` is::

    def reset_dq_bits(input,bits,extver=None,extname='dq')
    
These parameters now have to be described in a pair of configobj parameter files in order to define the parameters, their types and any validation that may need to be performed on the input values. 

Default Values for the Parameters
---------------------------------
The first file which needs to be defined provides the default values for each parameter.  Default values can be any string or numerical value, including "" or None.  

This task will simply need::

    _task_name_ = resetbits 
    input = "*flt.fits"
    bits = 4096
    extver = None
    extname = "dq"
 
The first line tells TEAL what task should be associated with this file. The default values for `extver` and `extname` simply match the defaults provided in the function signature. No default values were required for the other parameters, but these values were provided to support the most common usage of this task. 

This file needs to be saved with a filename extension of `.cfg` in a `pars/` subdirectory of the task's package. For `resetbits`, this file would be saved in the installation directory as the file::

    betadrizzle/lib/pars/resetbits.cfg
    
This file will then get installed in the directory `betadrizzle/pars/resetbits.cfg` with the instructions on how to set that up coming in the last step of this process.

Parameter Validation Rules
--------------------------
The type for the parameter values, along with the definition of any range of valid values, get defined in the second configobj file: the configobj specification (configspec) file or `cfgspc` file.  This file can also provide rules for how the GUI should respond to input values as well, turning the TEAL GUI into an active assistant for the user when editing large or complex sets of parameters. 

For this example, we have a very basic set of parameters to define without any advance logic required. TEAL provides validators for a wide range of parameter types, including:

  * `strings`: matches any string 
        Defined using `string_kw()` 
  * `integer`: matches any integer when a value is always required
        Defined using `integer_kw()`
  * `integer` or `None`: matches any integer or a value of None
        Defined using `integer_or_none_kw()` 
  * `float`: matches any floating point value, when a value is always required
        Defined using  `float_kw()` 
  * `float` or `None`: matches any floating point value or a value of None
        Defined using `float_or_none_kw()` 
  * `boolean`: matches boolean values - ``True`` or ``False``
        Defined using `boolean_kw()`
  * `option`: matches only those values provided in the list of valid options
        Defined using `option_kw()` command with the list of valid values as a parameter

There is also support for IP addresses as input parameters, and lists or tuples of any of these basic parameter types. Information on how to use those types, though, can be found within the `ConfigObj module`_ documentation.

With these available parameter types in mind, the parameters for the task can be defined in the configspec file. For the `resetbits` task, we would need::

    _task_name_ = string_kw(default="resetbits")
    input = string_kw(default="*flt.fits", comment="Input files (name, suffix, or @list)")
    bits = integer_kw(default=4096, comment="Bit value in array to be reset to 0")
    extver = integer_or_none_kw(default=None, comment="EXTVER for arrays to be reset")
    extname = string_kw(default="dq", comment= "EXTNAME for arrays to be reset")
    mode = string_kw(default="all")

Each of these parameter types includes a description of the parameter as the `comment` parameter, while default values can also be set using the `default` parameter value. This configspec file would then need to be saved alongside the .cfg file we just created as::
  
    betadrizzle/lib/pars/resetbits.cfgspc

.. note:: If you find that you need or want to add logic to have the GUI respond to various parameter inputs, this can always be added later by updating the parameter definitions in this file.  A more advanced example demonstrating how this can be done is provided in later sections. 


Step 2: TEAL Functions for the Task
===================================
TEAL requires that a couple of functions be defined within the task in order for the GUI to know how to get the help for the task and to run the task.  The functions that need to be defined are:

  * ``run(configobj=None)``
      This function serves as the hook to allow the GUI to run the task
  * ``getHelpAsString()``
      This function returns a long string which provides the help for the task

The sole input from TEAL will be a ConfigObj instance, a class which provides all the input parameters and their values after validation by the configobj validators.  This instance gets passed by TEAL to the tasks ``run()`` function and needs to be interpreted by that function in order to run the task.  

.. note:: The ``run()`` and ``getHelpAsString()`` functions, along with the task's primary user interface function, all need to be in module with the same name as the task as TEAL finds the task by importing the taskname. 

Defining the Help String
------------------------
The help information presented by the TEAL GUI comes from the ``getHelpAsString()`` function and gets displayed in a simple ASCII window.  The definition of this function can rely on help information included in the source code as docstrings or from an entirely separate file for tasks which have a large number of parameters or require long explanations to understand how to use the task.  The example from the `resetbits` task was simply::

    def getHelpAsString():
        helpString = 'resetbits Version '+__version__+__vdate__+'\n'
        helpString += __doc__+'\n'

        return helpString

This function simply relies on the module level docstring to describe how to use this task, since it is a simple enough task with only a small number of parameters. 

.. note:: The formatting for the docstrings or help files read in by this function can use the numpy documentation restructured text markup format to be compatible with Sphinx when automatically generating documentation on this task using Sphinx. The numpy extension results in simple enough formatting that works well in the TEAL Help window without requiring any translation. More information on this format can be found in the `Numpy Documentation`_ pages.

More complex tasks may require the documentation which would be too long to comfortably fit within docstrings in the code itself.  In those cases, separate files with extended discussions formatted using the numpy restructured text (reST) markup can be written and saved using the naming convention of ```<taskname>.help``` in the same directory as the module. The function can then simply use Python file operations to read it in as a list of strings which are concatenated together and passed along as the output.  The task `betadrizzle` uses separate files and can be used as an example of how this can be implemented. 


Defining How to Run the Task
----------------------------
The ConfigObj instance passed by TEAL into the module needs to be interpreted and used to run the application.  There are a couple of different models which can be used to define the interface between the ``run()`` function and the task's primary user interface function.  

  #. The ``run()`` function interprets the ConfigObj instance and calls the user interface  
     function. This works well for tasks which have a small number of parameters. 
     
  #. The ``run()`` function serves as the primary driver for the task and a separate 
     gets defined to provide a simpler interface for the user to call interactively. This
     works well for tasks which have a large number of parameters or sets of parameters
     defined in the ConfigObj interface. 
     
Our simple example for the task ``resetbits`` uses the first model, since it only has the 4 parameters as input. The ``run()`` function can simply be defined in this case as::

    def run(configobj=None):
        ''' Teal interface for running this code. '''

        reset_dq_bits(configobj['input'],configobj['bits'],
                      extver=configobj['extver'],extname=configobj['extname'])

    def reset_dq_bits(input,bits,extver=None,extname='dq'):

Interactive use of this function would use the function ``reset_dq_bits()``.  The TEAL interface would pass the parameter values in through the run function to parse out the parameters and send it to that same function as if it were run interactively. 


Step 3: Advertising TEAL-enabled Tasks 
======================================
Any task which has a TEAL interface implemented can be advertised to users of the package through the use of a ``teal`` function: ``teal.print_tasknames()``.  This function call can be added to the package's `__init__.py` module so that everytime the package gets imported, or reloaded, interactively, it will print out a message listing all the tasks which have TEAL GUI's available for use.  This listing will not be printed out when importing the package from another task.  The `__init__.py` module for the `betadrizzle` package has the following lines::

    # These lines allow TEAL to print out the names of TEAL-enabled tasks 
    # upon importing this package.
    from pytools import teal
    teal.print_tasknames(__name__, os.path.dirname(__file__))


Step 4: Setting Up Installation
===============================
The additional files which have been added to the package with the task now need to be installed alongside the module for the task.  Packages in the `STScI_Python` release get installed using Python's `distutils` mechanisms defined through the ``defsetup.py`` module. This file includes a dictionary for `setupargs` that describe the package and the files which need to be installed.  This needs to be updated to include all the new files as ``data_files`` by adding the following line to the ``setupargs`` dictionary definition::

  'data_files':  [(pkg+"/pars",['lib/pars/*']),( pkg, ['lib/*.help'])],
  
This will add the ConfigObj files in the `pars/` directory to the package while copying any ``.help`` files that were added to the same directory as the module. 


Step 5: Testing the GUI
=======================
Upon installing the new code, the TEAL interface will be available for the task.  There are a couple of ways of starting the GUI along with a way to grab the ConfigObj instance directly without starting up the GUI at all.

Running the GUI under PYRAF
---------------------------
The TEAL GUI can be started under PYRAF as if it were a standard IRAF task with the syntax::

    >>> import <package>
    >>> epar <taskname>
   
For example, our task ``resetbits`` was installed as part of the ``betadrizzle`` package, so we could start the GUI using::

    >>> import betadrizzle
    >>> epar resetbits
   
The fact that this task has a valid TEAL interface can be verified by insuring that the taskname gets printed out after the `import` statement.  

Running the GUI using Python 
----------------------------
Fundamentally, TEAL is a Python task that can be run interactively under any Python interpreter, not just PyRAF.  It can be called for our example task using the syntax::

    >>> from pytools import teal
    >>> cobj = teal.teal('resetbits')

Getting the ConfigObj Without Starting the GUI
----------------------------------------------
The function for starting the TEAL GUI, ``teal.teal()``, has a parameter to control whether or not to start the GUI at all.  The ConfigObj instance can be returned for the task without starting the GUI by using the `loadOnly` parameter. For our example task, we would use the command::

    >>> cobj = teal.teal('resetbits',loadOnly=True)
    
The output variable `cobj` can then be passed along or examined depending on what needs to be done at the time. 

---------------
Advanced Topics
---------------


.. _`ConfigObj module`: http://wiki.python.org/moin/ConfigObj
.. _`Numpy Documentation`: http://projects.scipy.org/numpy/wiki/CodingStyleGuidelines
