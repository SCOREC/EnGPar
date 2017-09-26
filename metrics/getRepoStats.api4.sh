#!/bin/bash -ex
jq=/users/cwsmith/software/jq-1.5/install/bin/jq
auth='Authorization: bearer <token>'
repourl=https://api.github.com/repos/SCOREC/EnGPar

stamp=`date +%d%m%Y`

#closed issues
curl https://api.github.com/graphql \
  -H "$auth" \
  -X POST \
  -o $stamp.issues \
  -d @- << EOF
{
  "query": "query {
     repository(owner:\"SCOREC\", name:\"EnGPar\") {
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

$jq '.data.repository.issues.edges[][]["updatedAt"]' $stamp.issues | xargs -L 1 date +"%d-%m-%Y" -d > engpar.issues
sort -t '-' -k 2 -k 1 engpar.issues > engpar.sorted.issues

#merged or closed pull requests 
curl https://api.github.com/graphql \
  -H "$auth" \
  -X POST \
  -o $stamp.pulls \
  -d @- << EOF
{
  "query": "query {
     repository(owner:\"SCOREC\", name:\"EnGPar\") {
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

$jq '.data.repository.pullRequests.edges[][]["updatedAt"]' $stamp.pulls | xargs -L 1 date +"%d-%m-%Y" -d > engpar.pulls
sort -t '-' -k 2 -k 1 engpar.pulls > engpar.sorted.pulls

#commits - this assumes that there are less than 100 commits per month
for year in 2015 2016 2017; do
  for i in {1..11}; do 
    echo $i
    printf -v mon "%02g" $i
    printf -v nextmon "%02g" $((i+1))
    echo 'mon' $mon 'nextmon' $nextmon
    log=$stamp.commits${year}.${mon}-${nextmon}
    curl https://api.github.com/graphql \
      -H "$auth" \
      -X POST \
      -o $log \
      -d @- << EOF
    {
      "query": "query {
         repository(owner:\"SCOREC\", name:\"EnGPar\") {
          ref(qualifiedName:\"master\") {
            name
            target {
              ... on Commit {
                history(first:100, since:\"${year}-${mon}-01T00:00:00Z\", until:\"${year}-${nextmon}-01T00:00:00Z\") {
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
      wc -l $log.name
      sort -u $log.name
    else
      echo "$log.name is empty"
      rm $log $log.name
    fi
  
  done
done
