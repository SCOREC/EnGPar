#!/bin/bash -ex
gnuplot \
  -e "inFile='core.views'" \
  -e "dateRangeStart='08092016'" \
  -e "dateRangeStop='25092017'" \
  -e "dataTitle='views'" \
  -e "outName='core.views.png'" \
  -e "yAxisLabel='unique views'" \
  -e "plotTitle='Unique Views of github.com/SCOREC/core'" \
  -e "idx=3" \
  -c plot.gplt

gnuplot \
  -e "inFile='core.clones'" \
  -e "dateRangeStart='08092016'" \
  -e "dateRangeStop='25092017'" \
  -e "dataTitle='clones'" \
  -e "outName='core.clones.png'" \
  -e "yAxisLabel='unique clones'" \
  -e "plotTitle='Unique Clones of github.com/SCOREC/core'" \
  -e "idx=3" \
  -c plot.gplt

gnuplot \
  -e "inFile='core.contributors'" \
  -e "dateRangeStart='08092016'" \
  -e "dateRangeStop='25092017'" \
  -e "dataTitle='contributors'" \
  -e "outName='core.contributors.png'" \
  -e "yAxisLabel='unique contributors'" \
  -e "plotTitle='Unique Contributors to github.com/SCOREC/core'" \
  -e "idx=2" \
  -c plot.gplt

gnuplot \
  -e "inFile='core.pulls'" \
  -e "dateRangeStart='08092016'" \
  -e "dateRangeStop='25092017'" \
  -e "dataTitle='pulls'" \
  -e "outName='core.pulls.png'" \
  -e "yAxisLabel='pull requests closed'" \
  -e "plotTitle='Pull Requests Closed github.com/SCOREC/core'" \
  -e "idx=2" \
  -c plot.gplt

gnuplot \
  -e "inFile='core.issues'" \
  -e "dateRangeStart='09092016'" \
  -e "dateRangeStop='25092017'" \
  -e "dataTitle='issues'" \
  -e "outName='core.issues.png'" \
  -e "yAxisLabel='issues closed'" \
  -e "plotTitle='Issues Closed github.com/SCOREC/core'" \
  -e "idx=2" \
  -c plot.gplt
