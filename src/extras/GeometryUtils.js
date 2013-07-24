/**
 * @author mrdoob / http://mrdoob.com/
 * @author alteredq / http://alteredqualia.com/
 */

var vIndex = ['a','b', 'c', 'd'];
var startIndex, endIndex, count, indexCount, compareIndex01, compareIndex02, compareEdge;
var edgeAlreadyProcessed = false;

function _sortAdjacentFaces(adjacentFacesArray){
    for (var i = 0; i < adjacentFacesArray.length; i++){
        //console.log("VERTEX NUMBER: " , i, "\n")
        _sort(adjacentFacesArray[i], i);
    }
}

function _sort(facesPerVertex, vertexIndex){
    startIndex = vertexIndex;
    var edges = [];

    for (var  i = 0; i < facesPerVertex.length; i++){

        facesPerVertex[edges.length].face instanceof THREE.Face4 ? count = 4 : count = 3;
        //find endIndex
        for (var j = 0; j < count; j++){
            if (facesPerVertex[edges.length].face[vIndex[j]] === startIndex){
                endIndex =   facesPerVertex[edges.length].face[vIndex[(j+1)%count]];
            }
        }
        //console.log("endIndex: ", endIndex);

        //if there aren't any edges in here yet, just push!
        if (edges.length < 1 ){
            edges.push( startIndex <= endIndex ?  {edge: [startIndex, endIndex] } : {edge: [endIndex, startIndex] } );
        }else{
            //else check if the edge is already in there! if it is, need to find new endIndex!
            for (var k = 0; k < edges.length; k++){
                if (! ((edges[k].edge[0]  === startIndex && edges[k].edge[1]=== endIndex) || (edges[k].edge[0]  === endIndex && edges[k].edge[1]=== startIndex )) )    {
                    edgeAlreadyProcessed = false;
                } else{
                    console.log("warning: edge already used");
                    //at this point, one would have to get a new endIndex!!!
                    edgeAlreadyProcessed = true;
                    break;
                }
            }
            if (!edgeAlreadyProcessed) {
                edges.push( startIndex <= endIndex ?  {edge: [startIndex, endIndex] } : {edge: [endIndex, startIndex] } );
            }
        }

        outerLoop:   for (j = edges.length; j < facesPerVertex.length; j++){
            //check if face is face3 or face4 and adjust loop-length appropriately
            facesPerVertex[j].face instanceof THREE.Face3 ? indexCount = 3: indexCount = 4;
            for ( k = 0; k < indexCount; k++){
                compareIndex01 =  facesPerVertex[j].face[vIndex[k]];
                compareIndex02 = facesPerVertex[j].face[vIndex[(k+1)%indexCount]];
                compareIndex01 <= compareIndex02 ? compareEdge = [compareIndex01, compareIndex02]  : compareEdge = [compareIndex02, compareIndex01];
                if (compareEdge[0] === edges[edges.length-1].edge[0] && compareEdge[1] === edges[edges.length-1].edge[1]){
                    //only swap if those 2 array elements do not follow each other in the array!
                    if ( j !== edges.length){
                        _swap(facesPerVertex, j, edges.length);
                        break outerLoop;
                    }
                }
            }
        }
    }
}

function _swap(array, index, startIndex ){
    // console.log("swap! ", index, startIndex)
    var temp = array[index];
    array[index] = array[startIndex];
    array[startIndex] = temp;
}

// calculates the average vector for one smoothing group by taking all non-normalized face normals (which are thereby area-weighted!) and just add them together.
//in the end, check if there are two smoothing groups that belong together, if so, add their vectors together and write them back on both averageVector attributes of the 2 smoothing groups
function _calculateVectorSumsForSmoothingGroups(connectedFaces, geometry){
    var smoothingGroupVectorSums = {};
    var smoothingGroup;

    for (var i = 0; i <connectedFaces.length; i++ ){
        //assign smoothing group ID by using only the first value of the array (in case it has more than 1 value)
        //in practice, as we where going "around the vertex", only the first/last face can have 2 smoothing groups
        smoothingGroup = connectedFaces[i].smoothingGroup[0];

        //if the groupID does not exist, use the newly calculated face normal as a start
        if ( !smoothingGroupVectorSums[smoothingGroup]){
            smoothingGroupVectorSums[smoothingGroup] = {averageVec: _calculateFaceNormal(connectedFaces[i].face, geometry)};

        }else{
            //check if groupID already used. if it is, instead of adding a new groupID entry, just add the un normalized face-normal
            smoothingGroupVectorSums[smoothingGroup].averageVec.add(_calculateFaceNormal(connectedFaces[i].face, geometry))
        }
    }

    //combine 2 smoothing groups if necessary
    if (connectedFaces[0].smoothingGroup.length === 2){

        var group1 = smoothingGroupVectorSums[connectedFaces[0].smoothingGroup[0]];
        var group2 = smoothingGroupVectorSums[connectedFaces[0].smoothingGroup[1]];

        var combinedVec = new THREE.Vector3().addVectors(group1.averageVec, group2.averageVec);

        //now both groups use the combined result
        group1.averageVec.copy(combinedVec);
        group2.averageVec.copy(combinedVec);
    }
    //console.log("sums: " ,smoothingGroupVectorSums)
    return smoothingGroupVectorSums;
}

function _calculateFaceNormal(face, geometry){
    //no normalization in this function as these are weighted normals
    var vA, vB, vC, vD;
    var cb = new THREE.Vector3(), ab = new THREE.Vector3(),
        db = new THREE.Vector3(), dc = new THREE.Vector3(), bc = new THREE.Vector3();

    if ( face instanceof THREE.Face3){
        vA = geometry.vertices[ face.a ];
        vB = geometry.vertices[ face.b ];
        vC = geometry.vertices[ face.c ];

        cb.subVectors( vC, vB );
        ab.subVectors( vA, vB );
        cb.cross( ab );
        return cb;
    }

    if (face instanceof  THREE.Face4){

        vA = geometry.vertices[ face.a ];
        vB = geometry.vertices[ face.b ];
        vC = geometry.vertices[ face.c ];
        vD = geometry.vertices[ face.d ];

        //triable abd
        db.subVectors( vD, vB );
        ab.subVectors( vA, vB );
        db.cross( ab );

        // triangle bcd
        dc.subVectors( vD, vC );
        bc.subVectors( vB, vC );
        dc.cross( bc );

        dc.add(db) ;
        return dc.multiplyScalar(0.5);
    }
}

THREE.GeometryUtils = {

	// Merge two geometries or geometry and geometry from object (using object's transform)

	merge: function ( geometry1, object2 /* mesh | geometry */, materialIndexOffset ) {

		var matrix, normalMatrix,
		vertexOffset = geometry1.vertices.length,
		uvPosition = geometry1.faceVertexUvs[ 0 ].length,
		geometry2 = object2 instanceof THREE.Mesh ? object2.geometry : object2,
		vertices1 = geometry1.vertices,
		vertices2 = geometry2.vertices,
		faces1 = geometry1.faces,
		faces2 = geometry2.faces,
		uvs1 = geometry1.faceVertexUvs[ 0 ],
		uvs2 = geometry2.faceVertexUvs[ 0 ];

		if ( materialIndexOffset === undefined ) materialIndexOffset = 0;

		if ( object2 instanceof THREE.Mesh ) {

			object2.matrixAutoUpdate && object2.updateMatrix();

			matrix = object2.matrix;

			normalMatrix = new THREE.Matrix3().getNormalMatrix( matrix );

		}

		// vertices

		for ( var i = 0, il = vertices2.length; i < il; i ++ ) {

			var vertex = vertices2[ i ];

			var vertexCopy = vertex.clone();

			if ( matrix ) vertexCopy.applyMatrix4( matrix );

			vertices1.push( vertexCopy );

		}

		// faces

		for ( i = 0, il = faces2.length; i < il; i ++ ) {

			var face = faces2[ i ], faceCopy, normal, color,
			faceVertexNormals = face.vertexNormals,
			faceVertexColors = face.vertexColors;

			if ( face instanceof THREE.Face3 ) {

				faceCopy = new THREE.Face3( face.a + vertexOffset, face.b + vertexOffset, face.c + vertexOffset );

			} else if ( face instanceof THREE.Face4 ) {

				faceCopy = new THREE.Face4( face.a + vertexOffset, face.b + vertexOffset, face.c + vertexOffset, face.d + vertexOffset );

			}

			faceCopy.normal.copy( face.normal );

			if ( normalMatrix ) {

				faceCopy.normal.applyMatrix3( normalMatrix ).normalize();

			}

			for ( var j = 0, jl = faceVertexNormals.length; j < jl; j ++ ) {

				normal = faceVertexNormals[ j ].clone();

				if ( normalMatrix ) {

					normal.applyMatrix3( normalMatrix ).normalize();

				}

				faceCopy.vertexNormals.push( normal );

			}

			faceCopy.color.copy( face.color );

			for ( var j = 0, jl = faceVertexColors.length; j < jl; j ++ ) {

				color = faceVertexColors[ j ];
				faceCopy.vertexColors.push( color.clone() );

			}

			faceCopy.materialIndex = face.materialIndex + materialIndexOffset;

			faceCopy.centroid.copy( face.centroid );

			if ( matrix ) {

				faceCopy.centroid.applyMatrix4( matrix );

			}

			faces1.push( faceCopy );

		}

		// uvs

		for ( i = 0, il = uvs2.length; i < il; i ++ ) {

			var uv = uvs2[ i ], uvCopy = [];

			for ( var j = 0, jl = uv.length; j < jl; j ++ ) {

				uvCopy.push( new THREE.Vector2( uv[ j ].x, uv[ j ].y ) );

			}

			uvs1.push( uvCopy );

		}

	},

	removeMaterials: function ( geometry, materialIndexArray ) {

		var materialIndexMap = {};

		for ( var i = 0, il = materialIndexArray.length; i < il; i ++ ) {

			materialIndexMap[ materialIndexArray[i] ] = true;

		}

		var face, newFaces = [];

		for ( var i = 0, il = geometry.faces.length; i < il; i ++ ) {

			face = geometry.faces[ i ];
			if ( ! ( face.materialIndex in materialIndexMap ) ) newFaces.push( face );

		}

		geometry.faces = newFaces;

	},

	// Get random point in triangle (via barycentric coordinates)
	// 	(uniform distribution)
	// 	http://www.cgafaq.info/wiki/Random_Point_In_Triangle

	randomPointInTriangle: function () {

		var vector = new THREE.Vector3();

		return function ( vectorA, vectorB, vectorC ) {

			var point = new THREE.Vector3();

			var a = THREE.Math.random16();
			var b = THREE.Math.random16();

			if ( ( a + b ) > 1 ) {

				a = 1 - a;
				b = 1 - b;

			}

			var c = 1 - a - b;

			point.copy( vectorA );
			point.multiplyScalar( a );

			vector.copy( vectorB );
			vector.multiplyScalar( b );

			point.add( vector );

			vector.copy( vectorC );
			vector.multiplyScalar( c );

			point.add( vector );

			return point;

		};

	}(),

	// Get random point in face (triangle / quad)
	// (uniform distribution)

	randomPointInFace: function ( face, geometry, useCachedAreas ) {

		var vA, vB, vC, vD;

		if ( face instanceof THREE.Face3 ) {

			vA = geometry.vertices[ face.a ];
			vB = geometry.vertices[ face.b ];
			vC = geometry.vertices[ face.c ];

			return THREE.GeometryUtils.randomPointInTriangle( vA, vB, vC );

		} else if ( face instanceof THREE.Face4 ) {

			vA = geometry.vertices[ face.a ];
			vB = geometry.vertices[ face.b ];
			vC = geometry.vertices[ face.c ];
			vD = geometry.vertices[ face.d ];

			var area1, area2;

			if ( useCachedAreas ) {

				if ( face._area1 && face._area2 ) {

					area1 = face._area1;
					area2 = face._area2;

				} else {

					area1 = THREE.GeometryUtils.triangleArea( vA, vB, vD );
					area2 = THREE.GeometryUtils.triangleArea( vB, vC, vD );

					face._area1 = area1;
					face._area2 = area2;

				}

			} else {

				area1 = THREE.GeometryUtils.triangleArea( vA, vB, vD ),
				area2 = THREE.GeometryUtils.triangleArea( vB, vC, vD );

			}

			var r = THREE.Math.random16() * ( area1 + area2 );

			if ( r < area1 ) {

				return THREE.GeometryUtils.randomPointInTriangle( vA, vB, vD );

			} else {

				return THREE.GeometryUtils.randomPointInTriangle( vB, vC, vD );

			}

		}

	},

	// Get uniformly distributed random points in mesh
	// 	- create array with cumulative sums of face areas
	//  - pick random number from 0 to total area
	//  - find corresponding place in area array by binary search
	//	- get random point in face

	randomPointsInGeometry: function ( geometry, n ) {

		var face, i,
			faces = geometry.faces,
			vertices = geometry.vertices,
			il = faces.length,
			totalArea = 0,
			cumulativeAreas = [],
			vA, vB, vC, vD;

		// precompute face areas

		for ( i = 0; i < il; i ++ ) {

			face = faces[ i ];

			if ( face instanceof THREE.Face3 ) {

				vA = vertices[ face.a ];
				vB = vertices[ face.b ];
				vC = vertices[ face.c ];

				face._area = THREE.GeometryUtils.triangleArea( vA, vB, vC );

			} else if ( face instanceof THREE.Face4 ) {

				vA = vertices[ face.a ];
				vB = vertices[ face.b ];
				vC = vertices[ face.c ];
				vD = vertices[ face.d ];

				face._area1 = THREE.GeometryUtils.triangleArea( vA, vB, vD );
				face._area2 = THREE.GeometryUtils.triangleArea( vB, vC, vD );

				face._area = face._area1 + face._area2;

			}

			totalArea += face._area;

			cumulativeAreas[ i ] = totalArea;

		}

		// binary search cumulative areas array

		function binarySearchIndices( value ) {

			function binarySearch( start, end ) {

				// return closest larger index
				// if exact number is not found

				if ( end < start )
					return start;

				var mid = start + Math.floor( ( end - start ) / 2 );

				if ( cumulativeAreas[ mid ] > value ) {

					return binarySearch( start, mid - 1 );

				} else if ( cumulativeAreas[ mid ] < value ) {

					return binarySearch( mid + 1, end );

				} else {

					return mid;

				}

			}

			var result = binarySearch( 0, cumulativeAreas.length - 1 )
			return result;

		}

		// pick random face weighted by face area

		var r, index,
			result = [];

		var stats = {};

		for ( i = 0; i < n; i ++ ) {

			r = THREE.Math.random16() * totalArea;

			index = binarySearchIndices( r );

			result[ i ] = THREE.GeometryUtils.randomPointInFace( faces[ index ], geometry, true );

			if ( ! stats[ index ] ) {

				stats[ index ] = 1;

			} else {

				stats[ index ] += 1;

			}

		}

		return result;

	},

	// Get triangle area (half of parallelogram)
	//	http://mathworld.wolfram.com/TriangleArea.html

	triangleArea: function () {

		var vector1 = new THREE.Vector3();
		var vector2 = new THREE.Vector3();

		return function ( vectorA, vectorB, vectorC ) {

			vector1.subVectors( vectorB, vectorA );
			vector2.subVectors( vectorC, vectorA );
			vector1.cross( vector2 );

			return 0.5 * vector1.length();

		};

	}(),

	// Center geometry so that 0,0,0 is in center of bounding box

	center: function ( geometry ) {

		geometry.computeBoundingBox();

		var bb = geometry.boundingBox;

		var offset = new THREE.Vector3();

		offset.addVectors( bb.min, bb.max );
		offset.multiplyScalar( -0.5 );

		geometry.applyMatrix( new THREE.Matrix4().makeTranslation( offset.x, offset.y, offset.z ) );
		geometry.computeBoundingBox();

		return offset;

	},

	triangulateQuads: function ( geometry ) {

		var i, il, j, jl;

		var faces = [];
		var faceUvs = [];
		var faceVertexUvs = [];

		for ( i = 0, il = geometry.faceUvs.length; i < il; i ++ ) {

			faceUvs[ i ] = [];

		}

		for ( i = 0, il = geometry.faceVertexUvs.length; i < il; i ++ ) {

			faceVertexUvs[ i ] = [];

		}

		for ( i = 0, il = geometry.faces.length; i < il; i ++ ) {

			var face = geometry.faces[ i ];

			if ( face instanceof THREE.Face4 ) {

				var a = face.a;
				var b = face.b;
				var c = face.c;
				var d = face.d;

				var triA = new THREE.Face3();
				var triB = new THREE.Face3();

				triA.color.copy( face.color );
				triB.color.copy( face.color );

				triA.materialIndex = face.materialIndex;
				triB.materialIndex = face.materialIndex;

				triA.a = a;
				triA.b = b;
				triA.c = d;

				triB.a = b;
				triB.b = c;
				triB.c = d;

				if ( face.vertexColors.length === 4 ) {

					triA.vertexColors[ 0 ] = face.vertexColors[ 0 ].clone();
					triA.vertexColors[ 1 ] = face.vertexColors[ 1 ].clone();
					triA.vertexColors[ 2 ] = face.vertexColors[ 3 ].clone();

					triB.vertexColors[ 0 ] = face.vertexColors[ 1 ].clone();
					triB.vertexColors[ 1 ] = face.vertexColors[ 2 ].clone();
					triB.vertexColors[ 2 ] = face.vertexColors[ 3 ].clone();

				}

				faces.push( triA, triB );

				for ( j = 0, jl = geometry.faceVertexUvs.length; j < jl; j ++ ) {

					if ( geometry.faceVertexUvs[ j ].length ) {

						var uvs = geometry.faceVertexUvs[ j ][ i ];

						var uvA = uvs[ 0 ];
						var uvB = uvs[ 1 ];
						var uvC = uvs[ 2 ];
						var uvD = uvs[ 3 ];

						var uvsTriA = [ uvA.clone(), uvB.clone(), uvD.clone() ];
						var uvsTriB = [ uvB.clone(), uvC.clone(), uvD.clone() ];

						faceVertexUvs[ j ].push( uvsTriA, uvsTriB );

					}

				}

				for ( j = 0, jl = geometry.faceUvs.length; j < jl; j ++ ) {

					if ( geometry.faceUvs[ j ].length ) {

						var faceUv = geometry.faceUvs[ j ][ i ];

						faceUvs[ j ].push( faceUv, faceUv );

					}

				}

			} else {

				faces.push( face );

				for ( j = 0, jl = geometry.faceUvs.length; j < jl; j ++ ) {

					faceUvs[ j ].push( geometry.faceUvs[ j ][ i ] );

				}

				for ( j = 0, jl = geometry.faceVertexUvs.length; j < jl; j ++ ) {

					faceVertexUvs[ j ].push( geometry.faceVertexUvs[ j ][ i ] );

				}

			}

		}

		geometry.faces = faces;
		geometry.faceUvs = faceUvs;
		geometry.faceVertexUvs = faceVertexUvs;

		geometry.computeCentroids();
		geometry.computeFaceNormals();
		geometry.computeVertexNormals();

		if ( geometry.hasTangents ) geometry.computeTangents();

	},

	setMaterialIndex: function ( geometry, index, startFace, endFace ){

		var faces = geometry.faces;
		var start = startFace || 0;
		var end = endFace || faces.length - 1;

		for ( var i = start; i <= end; i ++ ) {

			faces[i].materialIndex = index;

		}

    },
    calculateVertexNormals : function(geometry, angle){
        //reset vertex normals to zero-vectors
        for (var i = 0; i < geometry.faces.length; i++){
            if (geometry.faces[i] instanceof  THREE.Face3){
                geometry.faces[i].vertexNormals.length = 0;
                geometry.faces[i].vertexNormals.push( new THREE.Vector3() );
                geometry.faces[i].vertexNormals.push( new THREE.Vector3() ) ;
                geometry.faces[i].vertexNormals.push( new THREE.Vector3() );

            } else if (geometry.faces[i] instanceof  THREE.Face4 ){
                geometry.faces[i].vertexNormals.length = 0;
                geometry.faces[i].vertexNormals.push( new THREE.Vector3() );
                geometry.faces[i].vertexNormals.push( new THREE.Vector3() );
                geometry.faces[i].vertexNormals.push( new THREE.Vector3() );
                geometry.faces[i].vertexNormals.push( new THREE.Vector3() );
            }
        }

        //save face index per vertex index
        var adjacentNormals = [];
        var vN;

        for (var v = 0; v < geometry.vertices.length; v++){
            for (var f = 0; f < geometry.faces.length; f++){
                //this is needed for correct indexing of the given vertex with its vertexNormal.
                if (geometry.faces[f].a === v){
                    vN =  geometry.faces[f].vertexNormals[0];
                }else if(geometry.faces[f].b === v){
                    vN =  geometry.faces[f].vertexNormals[1];
                }else if(geometry.faces[f].c === v){
                    vN =  geometry.faces[f].vertexNormals[2];
                }else  if(geometry.faces[f].d === v){
                    vN =  geometry.faces[f].vertexNormals[3];
                }else{
                    vN = null;
                }
                if (vN !== null){
                    adjacentNormals[v] = adjacentNormals[v] || [];
                    adjacentNormals[v].push({face: geometry.faces[f], vertexNormal: vN, smoothingGroup:  []});
                }
            }
        }
        // sort faces in "adjacentNormals' because the face-objects are not sorted in an adjacent way, meaning that when going around 1 vertex by iterating through those faces
        // one does not know if face[i] and face[i+1] are adjacent in the array.
        _sortAdjacentFaces(adjacentNormals);

        //recalculate vertex normals
        var adjacentFaceNormal01 = new THREE.Vector3();
        var adjacentFaceNormal02 = new THREE.Vector3();
        var dotProduct;
        var angleBetweenFacesRad;
        var smoothing;
        var smoothingGroupID = 0;

        var normal, smoothingGroup;


        for (v = 0; v < geometry.vertices.length; v++){
            //for all faces the vertex is connected to
            for (i = 0; i < adjacentNormals[v].length; i++){
                //compare two adjacent faces (i) and (i+19 that are connected by the specified vertex v
                adjacentFaceNormal01.copy(adjacentNormals[v][i].face.normal);
                adjacentFaceNormal02.copy(adjacentNormals[v][(i+1)%adjacentNormals[v].length].face.normal);
                //calculate the dot product of the two face normals
                dotProduct = adjacentNormals[v][i].face.normal.dot( adjacentFaceNormal02 );

                //now calculate the angle between those 2 faces using the dot Product and the face normal length/ norm of the vector
                //result is in radian measure
                angleBetweenFacesRad = Math.acos( dotProduct / ( adjacentFaceNormal01.length() * adjacentFaceNormal02.length() ));

                if (THREE.Math.radToDeg(angleBetweenFacesRad) <= angle){
                    //console.log("angle:", _radToDeg(angleBetweenFacesRad)  + "    " ,i, (i+1)%adjacentNormals[v].length);
                    smoothing = true;

                    //if there are any smoothing groups, one will have to check if the array already contains the to be added groupID
                    if (adjacentNormals[v][i].smoothingGroup.length > 0){
                        for (var k = 0; k < adjacentNormals[v][i].smoothingGroup.length; k++){
                            if ( adjacentNormals[v][i].smoothingGroup[k] !== smoothingGroupID) {
                                adjacentNormals[v][i].smoothingGroup.push( smoothingGroupID);
                                break;
                            }
                        }
                    } else{
                        //else , smoothing group is still empty so just push the ID in
                        adjacentNormals[v][i].smoothingGroup.push( smoothingGroupID);
                    }
                    //same for the next (i+1) face in the array...
                    if (adjacentNormals[v][(i+1)%adjacentNormals[v].length].smoothingGroup.length > 0){
                        for ( k = 0; k < adjacentNormals[v][(i+1)%adjacentNormals[v].length].smoothingGroup.length; k++){
                            if ( adjacentNormals[v][(i+1)%adjacentNormals[v].length].smoothingGroup[k] !== smoothingGroupID) {
                                adjacentNormals[v][(i+1)%adjacentNormals[v].length].smoothingGroup.push(smoothingGroupID);
                                break;
                            }
                        }
                    } else{
                        adjacentNormals[v][(i+1)%adjacentNormals[v].length].smoothingGroup.push(smoothingGroupID);
                    }
                } else{
                    //if the angle is larger than specified, there is a hard edge between the 2 checked faces and therefore, increment smoothing group ID
                    //and use the existing face normal as the vertex normal
                    smoothingGroupID++;
                    adjacentNormals[v][i].vertexNormal.add(adjacentNormals[v][i].face.normal);
                    adjacentNormals[v][(i+1)%adjacentNormals[v].length].vertexNormal.add(adjacentNormals[v][(i+1)%adjacentNormals[v].length].face.normal);
                }
            }
            // console.log("counttest ", countTest);
            // console.log(adjacentNormals[v])


            // average vertexNormals based on smoothingGroupIDs
            //only loop if anything got smoothed at all
            if (smoothing){
                var groupIDNormals = _calculateVectorSumsForSmoothingGroups(adjacentNormals[v], geometry);

                for (k = 0; k < adjacentNormals[v].length; k++){
                    normal = adjacentNormals[v][k];
                    smoothingGroup = normal.smoothingGroup[0];
                    if (smoothingGroup !== undefined){
                        var groupNormal = groupIDNormals[smoothingGroup];
                        normal.vertexNormal.copy( groupNormal.averageVec );
                    }
                }
            }

            // reset values per vertex.
            // this also means that smoothingGroupIDs are local to their vertex with its connected faces.
            smoothing = false;
            smoothingGroupID = 0;
        }


        // if everything is calculated, last but not least normalize all vertex normals to have unit- vectors
        for (i = 0; i < geometry.faces.length; i++){
            for (var j = 0; j < geometry.faces[i].vertexNormals.length; j++){
                geometry.faces[i].vertexNormals[j].normalize();
            }
        }
    }

};
