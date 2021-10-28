	function mylen(string)

c	find true length of a string
	character string*(*)
	integer j, mylen

	do j=len(string),1,-1
		if (string(j:j).ne.' ') then
			mylen = j
			return
		end if
	end do
	mylen = j

	return
	end
