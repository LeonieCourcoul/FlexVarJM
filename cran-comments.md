## R CMD check results

0 errors | 0 warnings | 0 note

* This is a new release. I have fixed the notes about de new submission in replacing 'arxiv' by 'arXiv' in DESCRIPTION and ensuring that using my package do not change the user's options, par or working directory.
Concerning my examples, I have to estimate a model for each one and its takes time (a few minutes > 5 seconds), so it is impossible to have no notes without wrapped them in \donttest{}. However, when I run the check on my machine, the examples ran without notes. 
A big thank to the CRAN team for their check and feedback.
