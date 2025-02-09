Here we present a kind of coding/programming style policies to happily
coding together. First of all, as you may already noticed, we use
*C++* as our programming language. Why? well, that is because of its
widely continuous spreading along the scientific community and all its
nicely offered features.

The list of rules that you will find here are taken from our own
experience. We took what has worked for us when working in cooperative
projects, and set aside those silly things that only made cooperation
painful and no longer fun. We encourage you to talk with us about new
and excitement cooperation ideas you know. Hope you agree with us in
what we present here. Sit back, relax and happy reading.
  
* Comments

/"Comments are supposed to make your code easier to understand and
maintain-- not harder"/, [[http://blog.codinghorror.com/when-good-comments-go-bad][when comments go bad]].

There is nothing more frustrating than having to review someone's else
code (and sometimes our own code) and find out that the variables or
functions names are not representative of what they are used for. The
code is a complete mess, it has no identation at all, and the worse
thing is that there is not a single commented line to give us a clue
of what is going on. In this case we have two options,
- spend our whole day (or even more) to understand the code (and find
  out that the code is not doing what you are looking for), or
- throw it away and do your own implementation (I personally prefer
  this option).

Certainly, the above prevents code reusing which is one of the main
purposes of this project.

First of all, +try to+ make your code as *self-explanatory* as
possible. Suppose that you are telling a story to a group of people
and you want everyone that is hearing you do not get lost while you
are talking. If the syntaxis of the programming language does not
allow you to write your story clearly then use comments, do not be
afraid of using them.

Put your *comments as close as possible to the source code they are
referring*. A developer may not notice a comment that referrs to a
code line that he/she is modifying, thus the comments and the code
will be out of sync, [[http://blog.codinghorror.com/when-good-comments-go-bad][when commnets go bad]], [[http://blog.codinghorror.com/code-tells-you-how-comments-tell-you-why][code tells you how,
comments tell you why]].

Use *comments to clarify the selection of some algorithms* or to
mention things to remember during the execution of the code.

When writing *functions* that perform a large number of sub-tasks, use
the initial part of it to write a (brief!) *summary of the followed
strategy by the function* and include a list with the main steps of
the function. This helps the reader of your code to get an idea of
what to expect in the body of the function.

One final thing, always keep in mind that what may be transparently
obvious for you, may be completely obscure and complex to another
developer (or for yourself after returning from a non-programming
period), *use comments to ease the reading of your code, not to make
it harder.*

If you are looking for some good reasons for not documenting your code
then check these ones [[http://everything2.com/index.pl?node_id=1709851&displaytype=printable][(why programmers do not comment their code)]] and
see what better matches with you.

* Indentation
This is simple, *indent your code for easy reading*. I
use emacs as my preferred editor and have it configured to
automatically *indent with a single white-space* when pressing
TAB. You can copy and use my emacs configuration file (init.el) from
the */tools/emacs/* folder, or configure your favorite editor to
*indent with a single white-space*.

* Variables
We use variables to store values or data that are used frequently in
the body of a function or as part of a class.

** Variables types
Use =unsigned= instead of =int= in loops (=for= or =while=) that do
not require negative indexes.

** Variables names
Use *variable names that reflect the values or data stored in the
variable*. If you require to store the velocity you may use =v= or
=velocity= as the variable name.

Avoid using the well known variables names:

#+CAPTION: Bad variables names
#+NAME: tab:bad_variables_names
| =var1=            | =var2=           | =var3=             |
| =value=           | =variable=       | =my_variable=      |
| =this_variable=   | =other_variable= | =another_variable= |
| =delete_variable= | =here=           | =there=            |
| a                 | b                | c                  |

I know you (use?) have used those variable names, everyone does it,
but this is the time to forget about them and use self explained
variables names.

*Variables that stores the number of elements of something must use
the prefix* =n_=. For example, if you want to store the number
of processors in a variable then that variable must be named as
follow:
#+BEGIN_SRC c++
const unsigned n_processors = nprocessors();
#+END_SRC

Note that the function name does not use the '=_=', refer to the
functions names section for details.

** Const or non const?
Well, it happens that in C++, variables that do not pretend to change
their value along the entire exection of the program are declared with
a =const= before the variable type. Then why are they still called
variables?

Anyway, *use* =const= *on variables that are not intented to change
their value*. Remember that when using =const= you need to specify the
value of the variable at the time of its declaration.

* Functions
Functions are a great idea that let us split a complicated tasks in
small (or not that small) and easy to digest sub-tasks. We can
implement a complex task as a set of subtasks, each implementing a
basic idea that may be re-used in other complex tasks.

Think of a function as an independent task that may even call other
functions to perform its job.

** Functions types

** Functions names
When working in a small or individual project it is quite tempting to
use short name functions, first because no one else will use (or
review) our code, and second because of laziness. We pretend that this
library be (re-)used by a large community, thus function's names that
reflect the intention or the work performed by the function is a good
way to promote re-usability.

  - *Functions names MUST all be in lowercase*.
  * Use '=_=' to separate words in the function name.

** Split large funtions into sub-task

** Input and output arguments

A function may require some input data to work with, if that is the
case then you need to set it when calling the function. *Avoid using
global variables at all* to pass data to functions. Any function
should only know about the data that is receiving, if the function is
part of a class then the function should have access to the class
variables (including inherent data by the class).

*** Const or non-const 
Use =const= as much as you can, if you do not need (or do not know if
you need) to change the value of any variable inside a function then
use =const= after function arguments, example

#+BEGIN_SRC c++
  unsigned function_that_does_not_changes_values() const
#+END_SRC

otherwise do not use =const=
#+BEGIN_SRC c++
  unsigned function_that_does_changes_values()
#+END_SRC

Use =const= before the function name if the value that the function
returns is not expected (or if you dont know that it is expected) to
be modified by the function caller, example

#+BEGIN_SRC c++
  const unsigned function_whose_return_values_is_not_expected_to_change()
#+END_SRC

otherwise do not use =const=
#+BEGIN_SRC c++
  unsigned function_whose_return_values_is_expected_to_change()
#+END_SRC

*** Pass by copy or pass by reference
*Only pass arguments by copy when they are a single value*, such as an
integer or a double value. *Any other argument MUST be passed by
reference*. This is to avoid copying large vectors, matrices or
objects and thus run out of memory because of the many copies of the
same object in memory. If we do not really need a copy of every single
element in a vector, matrix or object then why should we make a copy
ot it?

Examples of passing arguments by reference here soon

Use \& when passing an argument by reference

* Classes

We use classes to represent entities that perform complex tasks, for
example, we use classes to implement linear algebra matrices. These
classes are in charge of providing storage, access and manipulation of
the matrices values. In order to identify abstract and concrete
classes we use the prefix =AC= for abstract classes and =CC= for
concrete classes. In general, abstract classes are used to define the
interfaces of the classes and the common data between
sub-classes. Concrete classes implement particular implementations of
the methods of the abstract class.

An abstract class for matrices is identified by the name =ACMatrix=,
and a concrete implementation of class representing matrices is
identified with the name =CCMatrix=.

** Member variables
** Member functions

DELETE DELETE

