import { BrowserRouter as Router, Switch, Route, Link as RouterLink } from 'react-router-dom';
import { Link, Button, Box } from '@material-ui/core';
import { useState, useEffect } from 'react';

function Menu() {
  const [fileName, setFilename] = useState('');

  useEffect(() => {
      setFilename(localStorage.getItem('fileName'));
  }, []);

  return (
    <Box
      display="flex"
      justifyContent="center"
      alignItems="center"
      minHeight="100vh"
    >
      <Button color="primary" component={RouterLink} to="/backendtest">
        BackendTest
      </Button>
      <Button color="secondary" component={RouterLink} to="/threetestphysics">
        ThreeTestPhysics
      </Button>
      <Button color="primary" component={RouterLink} to="/threetestgltf">
        ThreeTestGLTF
      </Button>
      <Button color="secondary" component={RouterLink} to="/threetest">
        ThreeTest {fileName}
      </Button>
    </Box>
  );
}

export default Menu;
