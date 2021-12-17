TODO
====

Checks
------

-   **(x)** Check that the dimensions/the content of the data given as
    inputs for the `predict` and the `simulate` methods are consistent
    with the content of the object used to predict or simulate from.

`summary` methods
-----------------

-   **(x)** Add a `summary` method for the class `"fgpm"`.

-   **(x)** Add a `summary` method for the class `"Xfgpm"`.

`plot` methods
--------------

-   Add a `plot` method for the class `"fgpm"`.

-   **(x)** Add a `plot` method for the class `"Xfgpm"`.

-   Add `autoplot` methods allowing to get \*\*ggplot2\* graphics.

Class `"Xfgpm"`
---------------

-   Store a copy of `sIn` and `sOut` as two new slots of the class. This
    will help at recreating any `fgpm` model visited from the data used
    in the fit. Maybe write a new function for that aim? How to specify
    the model to be re-created: by its number?

-   Improve the `summary` method. Should work correctly even when there
    are may inputs. Maybe add options to provide details on the scalar
    and the functional inputs separately in order to avoid having a
    table which is too wide.
