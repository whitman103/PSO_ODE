plot for[i=2:6] "trueOutData.txt" using (column(1)):(column(i)) with lines title "True 1",  "testOutData.txt" using (column(1)):(column(i))  with lines title "Test 1"


