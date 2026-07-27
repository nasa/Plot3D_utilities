import ReactDOM from "react-dom";
import { Canvas, useFrame, useLoader, useThree } from '@react-three/fiber';
import { OrbitControls, Stars, Text, Environment, GizmoHelper, GizmoViewport, useProgress, Html, Stats, PerspectiveCamera } from '@react-three/drei';
import { BoxGeometry, BufferGeometry, InstancedMesh, LineBasicMaterial, MeshNormalMaterial, SphereGeometry, ToneMapping, Line, ConvexHull } from "three";
import { ConvexGeometry } from 'three/examples/jsm/geometries/ConvexGeometry.js';
import { Suspense, useEffect, useState, useRef, useMemo, useCallback, useResource } from "react";
import * as THREE from 'three';
import usePromise from 'react-promise-suspense';
import { Line2 } from "three/examples/jsm/lines/Line2";
import DatGui, { DatBoolean, DatButton, DatColor, DatNumber, DatSelect, DatString, DatFolder, DatPresets } from "react-dat-gui";
import { useBeforeunload } from 'react-beforeunload';
import { useIdleTimer } from 'react-idle-timer';
import { saveAs } from 'file-saver';

import './Visualize.css';

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
      <boxGeometry args={[1*props.size, 1*props.size, 1*props.size]}/>
      <meshStandardMaterial color={props.color} opacity={props.visible ? props.opacity : 0.0} transparent />
    </instancedMesh>
  )
}

function dataToPoints (data) {
  let points = [];
  for (let i = 0; i < data.length; i++) {
    points.push(new THREE.Vector3(data[i][0], data[i][1], data[i][2]));
  }
  return points;
}


function lineFromPoints (points) {
  const geometry = new THREE.BufferGeometry().setFromPoints(points);
  const material = new THREE.LineBasicMaterial({ color: 0xffffff });
  return new THREE.Line(geometry, material);
}

function PointsToLines(props) {
  const ref = useRef();
  const points = [];
  for (let i = 0; i < props.coords.length; i++) {
    points.push(new THREE.Vector3(props.coords[i][0], props.coords[i][1], props.coords[i][2]));
  }
  const lineGeometry = new THREE.BufferGeometry().setFromPoints(points);
  return (
    <group>
      <line ref={ref} geometry={lineGeometry}>
        <lineBasicMaterial color={props.color} opacity={props.visible ? props.opacity : 0.0} transparent />
      </line>
    </group>
  );
}

function PointsToConvexGeometry(props) {
  const points = dataToPoints(props.coords);
  const geometry = new ConvexGeometry(points);
  const material = new MeshNormalMaterial( {color: 0x00ff00 } );
  return (
    <mesh geometry={geometry} material={material} />
  );
}

function Loader() {
  const { active, progress, errors, item, loaded, total } = useProgress()
  return <Html center>{progress} % loaded, waiting for API...</Html>
}

const fetchJson = input => fetch(input).then(res => res.json());


const CalcConnLabel = (calculating) => {
    if (calculating) {
        return("Calculating Connectivity...")
    }
    else {
        return("Calculate Connectivity")
    }
}

const CalcPeriodicityLabel = (calculating) => {
    if (calculating) {
        return("Calculating Periodicity...")
    }
    else {
        return("Calculate Periodicity")
    }
}

const CalcSplitblocksLabel = (calculating) => {
    if (calculating) {
        return("Calculating Split Blocks...")
    }
    else {
        return("Calculate Split Blocks")
    }
}

const randomColor = (colors) => {
  const randomColorChosen = colors[Math.floor(Math.random() * colors.length)];
  return randomColorChosen;
}

function Visualize() {
    const [blockData, setBlockData] = useState(null);
    const [connectivityData, setConnectivityData] = useState(null);
    const [periodicityData, setPeriodicityData] = useState(null);

    const [opts, setOpts] = useState({});
    const [blockOpts, setBlockOpts] = useState([]);
    const [connectivityOpts, setConnectivityOpts] = useState([]);
    const [periodicityOpts, setPeriodicityOpts] = useState([]);

    const [calculatingConn, setCalculatingConn] = useState(false);
    const [calculatingPeriodicity, setCalculatingPeriodicity] = useState(false);
    const [calculatingSplitblocks, setCalculatingSplitblocks] = useState(false);

    const [avgCentroidCoords, setAvgCentroidCoords] = useState(null);

    //const [target, setTarget] = useState(new THREE.Vector3(0, 0, 0));

    //const colors = ["silver", "gray", "maroon", "red", "purple", "fuchsia", "green", "lime", "olive", "yellow", "navy", "blue", "teal", "aqua"];

    const colors = ["red", "aqua", "white", "fuchsia"]

    const colors2 = ["green", "orange", "lime", "olive"]

    useEffect(() => {
      const requestOptions = {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({
              fileName: localStorage.getItem("fileName")
          })
      }

      fetch('/blocks', requestOptions)
          .then(response => response.json())
          .then(data => {
              // set data
              setBlockData(data);

              // set options for dat-gui
              const options = {
                  opacity: 1.0,
                  periodic_direction: "k",
                  rotation_axis: "x",
                  nblades: 55,
                  cellsPerBlock: 300000,
                  direction: "None",
              };
              Object.keys(data).forEach(k => {
                options[k] = true;
                options[k+"Color"] = colors[Number(k.slice(5, 6))];
                blockOpts.push(k);
              });
              setOpts(options);
              setBlockOpts(blockOpts);
          })
        
        fetch('/avgcentroid', requestOptions)
          .then(response => response.json())
          .then(data => {
              // set data
              setAvgCentroidCoords(new THREE.Vector3(data["cx_avg"], data["cz_avg"], data["cy_avg"])); // swapped coords of y and z because of how threejs works
          })
    }, []);

    const handleConnectivity = () => {
      const requestOptions = {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({
              fileName: localStorage.getItem("fileName")
          })
      };

      setCalculatingConn(true);

      fetch('/connectivities_grid', requestOptions)
          .then(response => response.json())
          .then(data => {
              setConnectivityData(data);
              const options = {};
              Object.keys(data).forEach(k => {
                  options[k] = true;
                  options[k+"Color"] = colors2[Number(k.slice(12, 13))];
                  // options[k+"Color"] = randomColor(colors);
                  connectivityOpts.push(k);
              });
              setOpts({...opts, ...options});
              setConnectivityOpts(connectivityOpts);
          })
          .then(() => setCalculatingConn(false));
    };

    const handlePeriodicity = () => {
      const requestOptions = {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body:JSON.stringify({
              fileName: localStorage.getItem("fileName"),
              periodicDirection: opts.periodic_direction,
              rotationAxis: opts.rotation_axis,
              nblades: opts.nblades
          })
      };

      setCalculatingPeriodicity(true);

      fetch('/periodicities', requestOptions)
          .then(response => response.json())
          .then(data => {
              setPeriodicityData(data);
              const options = {};
              Object.keys(data).forEach(k => {
                  options[k] = true;
                  options[k+"Color"] = randomColor(colors);
                  periodicityOpts.push(k);
              });
              setOpts({...opts, ...options});
              setPeriodicityOpts(periodicityOpts);
          })
          .then(() => setCalculatingPeriodicity(false));
    };

    const handleSplitblocks = () => {
      const requestOptions = {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body:JSON.stringify({
              fileName: localStorage.getItem("fileName"),
              cellsPerBlock: opts.cellsPerBlock,
              direction: opts.direction,
          })
        };

      setCalculatingSplitblocks(true);

      fetch('/splitblocks', requestOptions)
          .then(response => response.json())
          .then(data => {
              localStorage.setItem('fileName', data.fileName)
              setCalculatingSplitblocks(false);
                  
              })
      .then(() => window.location.reload()); 
    };

    const handleDownloadConnectivity = () => {
      const requestOptions = {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body:JSON.stringify({
              fileName: localStorage.getItem("fileName"),
              fileOption: "connectivity"
          })
        };

      fetch('/download', requestOptions)
        .then(response => response.blob())
        .then(data => saveAs(new Blob([data]), "connectivity.json"))
    };

    const handleDownloadPeriodicity = () => {
      const requestOptions = {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body:JSON.stringify({
              fileName: localStorage.getItem("fileName"),
              fileOption: "periodicity"
          })
        };

      fetch('/download', requestOptions)
        .then(response => response.blob())
        .then(data => saveAs(new Blob([data]), "periodicity.json"))
    };

    const handleDownloadSplitblocks = () => {
      const requestOptions = {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body:JSON.stringify({
              fileName: localStorage.getItem("fileName"),
              fileOption: "splitblocks"
          })
        };

      fetch('/download', requestOptions)
        .then(response => response.blob())
        .then(data => saveAs(new Blob([data]), "splitblocks.xyz"))
    };

    if ((blockData === null) && (avgCentroidCoords === null)) {
        return <h1>Loading, waiting for API...</h1>
    }

    return (
        <>
            <Canvas>
                <Stats />
                <ambientLight intensity={0.5} />
                <spotLight position={[0, 1, 0]} angle={0.3} />
                <Suspense fallback={<Loader />}>
                    <OrbitControls target={avgCentroidCoords} />
                    {blockOpts.map((opt) => (
                        <PointsToLines coords={blockData[opt]} visible={opts[opt]} opacity={opts.opacity} color={opts[opt+"Color"]} />
                    ))}
                    {connectivityOpts.map((opt) => (
                        <PointsToLines coords={connectivityData[opt]} visible={opts[opt]} opacity={opts.opacity} color={opts[opt+"Color"]} />
                    ))}
                    {periodicityOpts.map((opt) => (
                        <PointsToLines coords={periodicityData[opt]} visible={opts[opt]} opacity={opts.opacity} color={opts[opt+"Color"]} />
                    ))}
                </ Suspense>
                <GizmoHelper
                    alignment="bottom-right" // widget alignment within scene
                    margin={[80, 80]} // widget margins (X, Y)
                >
                    <GizmoViewport axisColors={['red', 'green', 'blue']} labelColor="black" />
                </GizmoHelper>
            </ Canvas>
            <DatGui data={opts} onUpdate={setOpts}>
                <DatNumber path="opacity" min={0.1} max={1.0} step={0.1} />
                <DatFolder title="blocks">
                    {blockOpts.map((opt) => (
                    <DatBoolean path={opt} />
                    ))}
                    {blockOpts.map((opt) => (
                    <DatSelect path={opt+"Color"} options={colors} />
                    ))}
                </ DatFolder>
                <DatButton label={CalcConnLabel(calculatingConn)} onClick={handleConnectivity} />
                <DatFolder title="Connectivities">
                    {connectivityOpts.map((opt) => (
                    <DatBoolean path={opt} />
                    ))}
                    {connectivityOpts.map((opt) => (
                    <DatSelect path={opt+"Color"} options={colors} />
                    ))}
                </ DatFolder>
                <DatFolder title="Periodicity Options">
                    <DatSelect path="periodic_direction" options={['i', 'j', 'k']} />
                    <DatSelect path="rotation_axis" options={['x', 'y', 'z']} />
                    <DatNumber path="nblades" min={1} step={1} />
                    <DatButton label={CalcPeriodicityLabel(calculatingPeriodicity)} onClick={handlePeriodicity} />
                </DatFolder>
                <DatFolder title="Periodicities">
                    {periodicityOpts.map((opt) => (
                    <DatBoolean path={opt} />
                    ))}
                    {periodicityOpts.map((opt) => (
                    <DatSelect path={opt+"Color"} options={colors} />
                    ))}
                </DatFolder>
                <DatFolder title="Split Blocks">
                    <DatNumber path="cellsPerBlock" min={1} step={1} />
                    <DatSelect path="direction" options={['i', 'j', 'k', "None"]} />
                    <DatButton label={CalcSplitblocksLabel(calculatingSplitblocks)} onClick={handleSplitblocks} />
                </DatFolder>
                <DatFolder title="Export to JSON">
                    <DatButton label="Export Connectivity" onClick={handleDownloadConnectivity} />
                    <DatButton label="Export Periodicity" onClick={handleDownloadPeriodicity} />
                    <DatButton label="Export Split Blocks" onClick={handleDownloadSplitblocks} />
                </DatFolder>
            </DatGui>
        </>
    );
}

export default Visualize;