#Interferometer Simulation Coding Convention

1. Function and Variable Naming
------------------------------------------------------------
* Use underscores to separate words:
	- Ex: `this_is_a_function`
	- Ex: `this_names_a_variable`

* Keep all letters lowercase

* Use verbose names:
	- The variables should describe their functionality
	- The function should describe its action
	    - Ex: a good name for a function that finds the perimeter of a square is:
	        `calculate_perimeter_square()`
	    - Ex: a good name for the perimeter of a square is:
	        `square_perimeter`

2. Comments
------------------------------------------------------------
* Block comments:
	- Keep in-line with current indentation
	- Use the following format:
	~~~~
		/* 
		 * First line
		 * Second line 
		 */
    ~~~~ 

* General comments:
	- If using a `//TODO`:
		- Give a description following the `TODO` on the same line
		- Provide a name on the `TODO` line indicating who should handle the task
	- Unless indicating a task to be completed, comments should be un-named
		- If the comment needs to be resolved and removed after resolution, attach a name to it
	- Make a tag name for `TODO` comments and indicate that tag at the start of the file central_simulator.c
		- Ex: Melanie Cornelius, tag: `Mcomment`

3. Pushing and Committing
------------------------------------------------------------
* Commit frequently, push only when information is needed by others or a major task has been completed
* Use verbose, descriptive commit messages
* Do not use versions, tags, or branches at this time
