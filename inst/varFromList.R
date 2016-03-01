mylist=(list(a=1,b=2,c="string1",d=list("r"=2,"z"="string2")))

for(i in 1:length(mylist)){
  ##first extract the object value
  tempobj=mylist[[i]]
  ##now create a new variable with the original name of the list item
  eval(parse(text=paste(names(mylist)[[i]],"= tempobj")))
}

# > print(a)
# [1] 1
# > print(b)
# [1] 2
# > print(c)
# [1] "string1"
# > print(d)
# $r
# [1] 2
# 
# $z
# [1] "string2"