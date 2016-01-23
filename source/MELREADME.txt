Hi Melanie,
This is Adam!

Before you look at the code in here, I actually want you to try running it. 
Change the Makefile so that the ROOT source points to where ROOT is on your system. Then you can just do 'make program' in the terminal to compile the program. 

Run it, and play along! PS: Gravity seems to be implemented correctly (phase shift from gravity is actually defined by one of Dr. Kaplan's articles), if I make the resolution high enough (~ 1000) I can see slight differences in the final interference pattern.

I know it's clunky as of right now (actually it's a lot better than before), but I was thinking that we can implement input files and have something to parse the input files. Input files are a big thing in the softwares used in my research, so I'm used to it. With as many parameters as we have to change, I think it's a necessity.

Now that you've maybe run the code a couple of times with different options:

Yes, you can look at the code now. The way I ended up making it work is probably really ugly; it's hard because a lot of the code is co-dependent. Most functions depend on other functions in other .cc files, However, it works! 

You'll notice plenty of global variables in the .c file. Some of them I feel justified about (const_e being the number e, pi, Plancks' constant, etc.), but some (like Gthick; thickness of the gratings) should end up being user-specified somehow. Same with cutoff, eta1 [open fraction of grating 1], eta2, period [of the gratings], etc.

I have a few tags for when you look through the code.
Search for 'varname' or 'Varname' when looking for bad variable names to change. 
Search for 'change' or 'efficient' or 'optimize' to see what I think could be changed. 

PS: It hit 2AM for me and I decided I needed some sleep.
