#!/bin/bash -e
jq=/path/to/jq/binary
auth='Authorization: token <your OAuth token>'
repourl=https://api.github.com/repos/SCOREC/EnGPar

stamp=`date +%d%m%Y`

#views
curl $repourl/traffic/views -H "$auth" -H "$media" > $stamp.views
count=$($jq -r '.count' $stamp.views)
unique=$($jq -r '.uniques' $stamp.views)
echo $stamp $count $unique >> core.views

#clones
curl $repourl/traffic/clones -H "$auth" -H "$media" > $stamp.clones
count=$($jq -r '.count' $stamp.clones)
unique=$($jq -r '.uniques' $stamp.clones)
echo $stamp $count $unique >> core.clones

#send a stats request to ensure stats have been compiled
curl $repourl/stats/commit_activity -H "$auth" > /dev/null
sleep 60 # give github time to compile stats

#weekly commits - collect previous weeks stats on a monday
dayofweek=$(date +%u)
if [ $dayofweek = 1 ]; then # day 1 is monday
  curl $repourl/stats/commit_activity -H "$auth" > $stamp.commits
  count=$($jq -r ".[50].total" $stamp.commits) # week 50 is the previous week
  echo $stamp $count >> core.commits
fi

#unique committers
curl $repourl/stats/contributors -H "$auth" > $stamp.contributors
count=$($jq -r ".[].total" $stamp.contributors | wc -l)
echo $stamp $count >> core.contributors

#pull requests
curl $repourl/pulls?state=closed -H "$auth" > $stamp.pulls
count=$($jq -r ".[].id" $stamp.pulls | wc -l)
echo $stamp $count >> core.pulls

#closed issues
curl $repourl/issues?state=closed -H "$auth" > $stamp.issues
count=$($jq -r ".[].title" $stamp.issues | wc -l)
echo $stamp $count >> core.issues

rm $stamp.*
