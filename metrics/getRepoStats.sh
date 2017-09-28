#!/bin/bash -x

cd /users/cwsmith/Projects/engpar/

jq=/users/cwsmith/software/jq-1.5/install/bin/jq
auth='Authorization: bearer <token>'
media='Accept:   application/vnd.github.spiderman-preview'

stamp=`date +%d%m%Y`

#weekly commits - this assumes that there are less than 100 commits per week
getCommitsAndDevs() {
  local repo=$1
  local branch=$2
  local dayofweek=$(date +%u)
  if [ $dayofweek = 4 ]; then
    local startdate=$(date +%Y-%m-%d -d -1week)
    local enddate=$(date +%Y-%m-%d)
    echo "startdate \'$startdate\' enddate \'$enddate\'"
    local log=$stamp.commits.${startdate}_${enddate}
    curl https://api.github.com/graphql \
      -H "$auth" \
      -X POST \
      -o $log \
      -d @- << EOF
{
  "query": "query {
     repository(owner:\"SCOREC\", name:\"${repo}\") {
      ref(qualifiedName:\"${branch}\") {
        name
        target {
          ... on Commit {
            history(first:100, since:\"${startdate}T00:00:00Z\", until:\"${enddate}T00:00:00Z\") {
              edges {
                node {
                  oid
                  message
                  author {
                    name
                    email
                    date
                  }
                }
              }
            }
          }
        }
      }
    }
  }"
}
EOF
    
    $jq '.data.repository.ref.target.history.edges[][]["author"]["name"]' $log > $log.name
    if [ -s $log.name ]; then
      local count=$(wc -l < $log.name)
      echo $stamp $count >> $repo.commits
      local devs=$(sort -u $log.name | wc -l)
      echo $stamp $devs >> $repo.devs
    fi
  fi
}

#closed issues - assumes less than 100 closed issues in the last week
getClosedIssues() {
  local repo=$1
  local dayofweek=$(date +%u)
  if [ $dayofweek = 4 ]; then
    local lastweek=$(date +%Y-%m-%d -d -1week)
    local today=$(date +%Y-%m-%d)
    local log=$stamp.issues.${lastweek}_${today}
    curl https://api.github.com/graphql \
      -H "$auth" \
      -X POST \
      -o $log \
      -d @- << EOF
{
  "query": "query {
     repository(owner:\"SCOREC\", name:\"${repo}\") {
       issues(last:100, states:CLOSED) {
         edges {
           node {
             title
             updatedAt
           }
         }
       }
     }
  }"
}
EOF

    if [ -s $log ]; then
      $jq ".data.repository.issues.edges[][] | select ( .updatedAt > \"${lastweek}T00:00:00Z\" )" $log > $log.lastweek
      count=$(grep "updatedAt" $log.lastweek)
      echo $stamp $count >> $repo.issues
    fi
  fi
}

#merged or closed pull requests - assumes less than 100 in the last week
getMergedOrClosedPullRequests() {
  local repo=$1
  local dayofweek=$(date +%u)
  if [ $dayofweek = 4 ]; then
    local lastweek=$(date +%Y-%m-%d -d -1week)
    local today=$(date +%Y-%m-%d)
    local log=$stamp.pulls.${lastweek}_${today}
    curl https://api.github.com/graphql \
      -H "$auth" \
      -X POST \
      -o $log \
      -d @- << EOF
{
  "query": "query {
     repository(owner:\"SCOREC\", name:\"${repo}\") {
       pullRequests(last:100, states:[CLOSED, MERGED]) {
         edges {
           node {
             title
             state
             updatedAt
             mergedAt
           }
         }
       }
     }
  }"
}
EOF

    if [ -s $log ]; then
      $jq ".data.repository.pullRequests.edges[][] | select ( .updatedAt > \"${lastweek}T00:00:00Z\" )" $log > $log.lastweek
      count=$(grep "updatedAt" $log.lastweek)
      echo $stamp $count >> $repo.pulls
    fi
  fi
}

##### core #####

repourl=https://api.github.com/repos/SCOREC/core
repo=core
branch=develop

#views
curl $repourl/traffic/views -H "$auth" -H "$media" > $stamp.views
count=$($jq -r '.count' $stamp.views)
unique=$($jq -r '.uniques' $stamp.views)
echo $stamp $count $unique >> $repo.views

#clones
curl $repourl/traffic/clones -H "$auth" -H "$media" > $stamp.clones
count=$($jq -r '.count' $stamp.clones)
unique=$($jq -r '.uniques' $stamp.clones)
echo $stamp $count $unique >> $repo.clones

#commits and devs
getCommitsAndDevs "$repo" "develop"

#closed issues
getClosedIssues "$repo"

#merged or closed pull requests
getMergedOrClosedPullRequests "$repo"

##### EnGPar #####

repourl=https://api.github.com/repos/SCOREC/EnGPar
repo=EnGPar
branch=master

#views
curl $repourl/traffic/views -H "$auth" -H "$media" > $stamp.views
count=$($jq -r '.count' $stamp.views)
unique=$($jq -r '.uniques' $stamp.views)
echo $stamp $count $unique >> engpar.views

#clones
curl $repourl/traffic/clones -H "$auth" -H "$media" > $stamp.clones
count=$($jq -r '.count' $stamp.clones)
unique=$($jq -r '.uniques' $stamp.clones)
echo $stamp $count $unique >> engpar.clones

#commits and devs
getCommitsAndDevs "$repo" "$branch"

#closed issues
getClosedIssues "$repo"

#merged or closed pull requests
getMergedOrClosedPullRequests "$repo"

#rm $stamp.*
