#!/bin/bash
export DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
if ! [ "${GITHUB_ACTIONS}" == "true" ]; then
	>&2 echo "ERROR: Not running on Travis (GITHUB_ACTIONS='${GITHUB_ACTIONS}')."
	exit 1
fi

# Check that the documentation has indeed been built
if ! [ -d "${DIR}/../doc/html" ]; then
	>&2 echo "ERROR: documentation not built (doc/html missing). Aborting."
	tree "${DIR}/.." # for debugging
	exit 1
fi

# Write the SSH deploy key to a file.
if [ -z ${GITHUB_DEPLOY_KEY+x} ]; then
	>&2 echo "ERROR: \$GITHUB_DEPLOY variable is unset."
	exit 1
fi
mkdir -p ~/.ssh
base64 -d > ~/.ssh/github_deploy_key <<-EOF
	${GITHUB_DEPLOY_KEY}
EOF
chmod og-rwx ~/.ssh/github_deploy_key # the SSH key needs conservative permissions

# We'll set up a custom SSH config file and tell Git to use it with the GIT_SSH_COMMAND
# environment variable.
cat - >> ~/.ssh/grasp_deploy_config <<-EOF
Host github.com
	StrictHostKeyChecking no
	HostName github.com
	IdentityFile ~/.ssh/github_deploy_key
EOF
export GIT_SSH_COMMAND="ssh -F ~/.ssh/grasp_deploy_config"

# Pushing the documentation to the gh-pages branch. First, we'll create a temporary
# directory and clone the repository again, and cd there:
gh_pages_dir=`mktemp -d`
echo "INFO: Cloning gh-pages into ${gh_pages_dir}"
git clone -b gh-pages git@github.com:${GITHUB_REPOSITORY} "${gh_pages_dir}" || exit
cd "${gh_pages_dir}" || exit
# Remove all the old files (but not stuff under .git)
find . -type f -not -regex "\./\.git/.*" -delete
# Copy over the new docs:
rsync -a ${DIR}/../doc/html/ . || exit
# Git add, commit and push the updated documentation:
git add -A || exit
git commit -m "Deploy documentation" || {
	echo "git commit returned a non-zero exit code ($?). This is expected if there are no changes to commit."
	exit 0
}
git push --set-upstream origin gh-pages || exit
echo "Documentation deployment complete."
