## A tip on the shebang line

If the '<avl>' does not run, it may be caused by missing a shebang line at the beginning.
The issue can be solved by adding a shebang line to the beginning of the file as described in https://unix.stackexchange.com/questions/491419/fish-shell-exec-format-error.

In our case, if the error does occur, try to invoke python, and then type

```
<with open('avl','w') as f: f.write('#!/bin/sh\nexit 0')>
```
