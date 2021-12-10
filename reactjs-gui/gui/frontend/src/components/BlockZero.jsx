import React from "react";
import { fetchBlockZeroURL } from "./Endpoints";

import { useRef, useEffect } from "react";
import * as THREE from 'three';

const resource = fetchBlockZeroURL()

// Use instanceMesh to greatly improve performance when drawing many boxes of same geometry and material
function InstancedBoxes(props) { // { temp = new THREE.Object3D() }) broken because data doesnt com in fast enough
    // NEED TO LEARN HOW TO USE PROMISES?
    const ref = useRef();
    useEffect(() => {
      // Set positions
      console.log(props.coords.length)
      for (let i = 0; i < props.coords.length; i++) {
        props.temp.position.set(props.coords[i][0], props.coords[i][1], props.coords[i][2])
        props.temp.updateMatrix()
        ref.current.setMatrixAt(i, props.temp.matrix)
      }
      // Update the instance
      ref.current.instanceMatrix.needsUpdate = true
    }, [])
    return (
      <instancedMesh ref={ref} args={[null, null, props.coords.length]}>
        <boxGeometry args={[0.05, 0.05, 0.05]}/>
        <meshStandardMaterial color="red" />
      </instancedMesh>
    )
}

const BlockZero = () => {
    const coords = resource.read()

    return (
        <InstancedBoxes coords={coords} color="red" temp={ new THREE.Object3D() } />
    )
}

export default BlockZero