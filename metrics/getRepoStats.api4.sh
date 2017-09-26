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

$jq '.data.repository.issues.edges[][]["updatedAt"]' 26092017.issues | xargs -L 1 date +"%d%m%Y" -d > engpar.issues
sort -t '-' -k 2 -k 1 engpar.issues > engpar.sorted.issues

#rm $stamp.*
