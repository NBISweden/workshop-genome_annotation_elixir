# workshop-common, Elixir mode.
This version of the common coursepage layout uses ELIXIR instead of the NBIS as the main logo.

### How to use it
If you want to create a new course page, first follow the common guidelines
given [here](https://nbisweden.github.io/workshop-howto/labs/new_course).


Then tell git to use this branch:
```
cd </path/to/your/workshop-tba>
git config -f .gitmodules submodule.workshop-common.branch feature/elixir_mode
git submodule update --remote
```
Add and commit this change:
```
git add workshop-common .gitmodules
git commit -m "Use elixir layout"
git push
```
After a few minutes, your github page should be updated to reflect the changes. 
