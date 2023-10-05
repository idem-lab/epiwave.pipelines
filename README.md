# GPreff

Repository for the new $R_eff$ model using Gaussian Processes with `greta`. This model is centered on estimation of the infection timeseries, from which $R_eff$ can be calculated when we add generation interval information.

This workflow has separate components for PCR and RAT data, and estimates separately by jurisdiction, in this case Australian state.

## Instructions for contributing to this repository
### Setting up

1. Fork this repo into personal GitHub account
2. Navigate to your fork [username]/GPreff and clone a local version to your machine
3. Open new project from version control in local RStudio
4. Configure your remote repository

  * list remotes
   ```
   git remote -v
   > origin  git@github.com:[YOURUSERNAME]/GPreff.git (fetch)
   > origin  git@github.com:[YOURUSERNAME]/GPreff.git (push)
   ```

  * add upstream
   ```
   git remote add upstream git@github.com:infectious-disease-ecology-modelling/GPreff.git
   ```

  * verify by checking remotes again
  ```
  git remote -v
   > origin  git@github.com:[YOURUSERNAME]/GPreff.git (fetch)
   > origin  git@github.com:[YOURUSERNAME]/GPreff.git (push)
   > upstream        git@github.com:infectious-disease-ecology-modelling/GPreff.git (fetch)
   > upstream        git@github.com:infectious-disease-ecology-modelling/GPreff.git (push)
  ```

### Making changes

1. Before beginning new work, make sure your fork is up to date with the remote

   ```
   # if you are not on main branch
   git checkout main

   # fetch from upstream repo
   git fetch upstream main

   # then merge
   git merge upstream/main
   ```

2. Make a new branch for your changes using:

   ```
   git checkout -b ＜dev-branch＞
   ```
   Try and use a descriptive name for the new branch, for example 'updateCAR'

3. When you've made your changes, commit all changes to your branch, then repeat step 1 if necessary (if there have been changes to upstream)

4. Merge your new branch into your main, using `--no-ff` to indicate 'no fast forward' because we want a commit recording the merge

   ```
   git checkout main
   git merge --no-ff <dev-branch>
   ```
  
5. Add commit message, then push to origin and delete dev branch

   ```
   git push origin main
   git branch --delete <dev-branch>
   ```

5. Create pull request on GitHub


### Resources
* For more on working with branches, see [here](https://www.freecodecamp.org/news/how-to-work-with-branches-in-git/).
* For more on forks and working with upstreams, see [here](https://www.atlassian.com/git/tutorials/git-forks-and-upstreams).
* For more on merging and rebasing, see [here](https://www.atlassian.com/git/tutorials/merging-vs-rebasing).

