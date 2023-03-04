Point cloud: global scope

unsigned int Function expand_d-1_simplex (input: pair(d-1_face, point))
- Get circumcenter(`p1`) of `d-1_face`
- Get circumcenter(`p2`) of `d-1_face`+`point`
- Create halfspaces using p1 and p2 and partition the entire pointcloud using normal vector to plane(`p1`&`p2`) and vector(`data_point`-`p1`) // Use a.b
- if `opposite_half_space`==NULL
  -return -1 //Triangulation point can only exist in half space opposite to the context_point
- Store circumradius for all `d-1_face`+`data_point`(opposte_half_space) in vector.
- Find `minimum radius point`
- Find circumcenter with minimum radius point. Find points inside minimum radius circle.
- Find point with maximum circumradius point(`triangulation_point`) inside the minimum radius circle.
- return `triangulation_point`


Function main()
- Declare a set of set called `complete_delaunay`
- Initialize `current_delaunay` with an initial heuristic //Later moved to a seperate function(found by serial insertion of points into single intial triangle)
- Initialize a map(set,unsigned) of d-1 simplices and context point (combination was previous triangulation) of `current_context` called // outermost skin of the current complete delaunay
- Declare map(set,unsigned) of d-1 simplices and context point (combination was previous triangulation) of `next_context`

- while (`current_context` is not null)
    - `next_context` << NULL
    - For each pair(key,value) in `current_context` // Run in Parallel(Map-reduce)
      - Run `expand_d-1_simplex` on the pairs  //Map task
      - If a point is returned that is not -1
        - Add all combinations of `d-1_simplex` and the returned point to `next_context` //Reduce task
        //Do not include collided faces(>1 returns for a face) && faces from `current_context` again
    - Push triangulations from `next_context` to `complete_delaunay`
    - `current_context` << `next_context`