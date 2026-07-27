import React, { useRef, useFrame } from 'react'
import { useGLTF, PerspectiveCamera } from '@react-three/drei'

export default function Model(props) {
  const group = useRef()
  const { nodes, materials } = useGLTF('/finalmesh_points_block1_gltf.gltf')
  return (
    <group ref={group} {...props} dispose={null}>
      <PerspectiveCamera
        makeDefault={false}
        far={88.8}
        near={60.55}
        fov={30}
        position={[0.58, 8.29, -55.59]}
        rotation={[3.04, -0.15, 2.08]}
      />
      <mesh
        castShadow
        receiveShadow
        geometry={nodes.mesh0.geometry}
        material={nodes.mesh0.material}
      />
      <mesh
        castShadow
        receiveShadow
        geometry={nodes.mesh1.geometry}
        material={nodes.mesh1.material}
      />
    </group>
  )
}

useGLTF.preload('/finalmesh_points_block1_gltf.gltf')

