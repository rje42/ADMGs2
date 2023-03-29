##' Get m-connected vertices
##' 
##' @param graph an object of class \code{mixedgraph}
##' @param X set of vertices to check for m-connection
##' @param Z vertices to condition on
##' 
##' @details Perform a breadth-first search to find m-connected vertices given 
##' a specified set.
##' 
##' Based on the C++ function 'dConnected' from the package \code{dagitty}.
##' 
##' @export
m_connected <- function(graph, X, Z) {  # Y, AnZ
	# graph2 <- graph	
	
	## check only allowed edges are included
	allow_edges <- match(c("undirected", "directed", "bidirected"
                        # "partially directed", "partially undirected"
	                       ), 
	                     edgeTypes()$type, nomatch = 0L)
	
	if (nedge(graph, edges = edgeTypes()$type[-allow_edges]) > 0) {
	  stop(paste0("Cannot have edges of type ", paste(edgeTypes()$type[-allow_edges], collapse=", ")))
	}
	
	# if (nedge(graph, edges = "partially directed") > 0) {
	#   graph$edges$directed <- collapse(graph$edges[c("directed", "partially directed")])
	# }
	# if (nedge(graph, edges = "partially undirected") > 0) {
	#   graph$edges$directed <- collapse(graph$edges[c("directed", "partially undirected")], dir=c(1,-1))
	#   collapse(graphCr("1 -> 2 o- 1")$edges[c("directed", "partially undirected")], dir=c(1, 1))
	# }


	# if (missing(anZ)) anZ <- anc(graph, Z)
	## get ancestors of Z
	anZ <- anc(graph, Z)
	
	frwd_q <- frwd_v <- bkwd_v <- integer(0)
	bkwd_q <- X
		# var forward_queue = []
		# var backward_queue = []
		# var forward_visited ={}
		# var backward_visited = {}
		# var i, Y_ids = {}, Z_ids = {}, AnZ_ids = {}, v, vv

		# for( i = 0 ; i < X.length ; i ++ ){
		# 	backward_queue.push( X[i] )
		# }
		# for( i = 0 ; i < Y.length ; i ++ ){
		# 	Y_ids[Y[i].id] = 1
		# }
		# for( i = 0 ; i < AnZ.length ; i ++ ){
		# 	AnZ_ids[AnZ[i].id] = 1
		# }
		# for( i = 0 ; i < Z.length ; i ++ ){
		# 	Z_ids[Z[i].id] = 1
		# }

		while ( length(frwd_q) + length(bkwd_q) > 0 ) {
			if ( length(frwd_q) > 0 ) {
				v <- frwd_q[1]
				frwd_q <- frwd_q[-1]
				frwd_v <- c(frwd_v, v)
				
				# if( Y_ids[v.id] ) return true
				if( v %in% anZ ) {
				  ## add parents to backwards queue
					pas_v <- pa(graph, v)
					bkwd_q <- c(bkwd_q, setdiff(pas_v, bkwd_v))
					# for( i = 0 ; i < vv.length ; i ++ ){
					# 	if( !backward_visited[vv[i].id] ){
					# 		backward_queue.push( vv[i] )
					# 	}
					# }
					
					## add siblings to forwards queue
					sib_v <- sib(graph, v)
					frwd_q <- union(frwd_q, setdiff(sib_v, frwd_v))
					# for( i = 0 ; i < vv.length ; i ++ ){
					# 	if( !forward_visited[vv[i].id] ){
					# 		forward_queue.push( vv[i] )
					# 	}
					# }
				} 
				if( !(v %in% Z) ){
				  ## add children and neighbours to forward queue
				  chn_v <- c(ch(graph, v), nb(graph, v))
				  frwd_q <- union(frwd_q, chn_v)
				  
					# vv = _.union( v.getChildren(), v.getNeighbours() )
					# for( i = 0 ; i < vv.length ; i ++ ){
					# 	if( !forward_visited[vv[i].id] ){
					# 		forward_queue.push( vv[i] )
					# 	}
					# }
				}
			}
			if( length(bkwd_q) > 0 ) {
				v <- bkwd_q[1]
				bkwd_q <- bkwd_q[-1]
				bkwd_v <- c(bkwd_v, v)
				
				# if( Y_ids[v.id] ) return true
				if(v %in% Z) next
				
				## add children and siblings to forward queue
				chsb_v <- c(ch(graph, v), sib(graph, v))
				frwd_q <- c(frwd_q, setdiff(chsb_v, frwd_v))
				
				# vv = _.union( v.getChildren(), v.getSpouses() )
				# for( i = 0 ; i < vv.length ; i ++ ){
				# 	if( !forward_visited[vv[i].id] ){
				# 		forward_queue.push( vv[i] )
				# 	}
				# }
				## add parents and neighbours to backward queue
				panb_v <- c(pa(graph, v), nb(graph, v))
				bkwd_q <- c(bkwd_q, setdiff(panb_v, bkwd_v))
				# 
				# vv = _.union( v.getParents(), v.getNeighbours() )
				# for( i = 0 ; i < vv.length ; i ++ ){
				# 	if( !backward_visited[vv[i].id] ){
				# 		backward_queue.push( vv[i] )
				# 	}
				# }
			}
		}
		# if( Y.length > 0 ){
		# 	return false
		# } else {
		# 	return go.getVertex(
		# 		_.union( Object.keys( forward_visited ), Object.keys( backward_visited ) ) 
		# 	)
		# }
	out <- union(frwd_v, bkwd_v)
	return(out)	
}