import './App.css';
import { BrowserRouter as Router, Switch, Route, Link as RouterLink } from 'react-router-dom';
import Menu from './pages/Menu';
import Upload from './pages/Upload';
import Visualize from './pages/Visualize';

import { useMemo } from 'react';
import useMediaQuery from '@mui/material/useMediaQuery';
import { createTheme, ThemeProvider } from '@mui/material/styles';
import CssBaseline from '@mui/material/CssBaseline';

import { useIdleTimer } from 'react-idle-timer';

// TODO: update to react router v6
// https://reactrouter.com/docs/en/v6/upgrading/v5#upgrade-to-react-v168
function App() {
  const prefersDarkMode = useMediaQuery('(prefers-color-scheme: dark)');

  const theme = useMemo(
    () =>
      createTheme({
        palette: {
          mode: prefersDarkMode ? 'dark' : 'light',
        },
      }),
    [prefersDarkMode],
  );

  // Session timeout after 2 hours
  const IDLE_MINUTES = 120;

  const handleOnIdle = event => {
    // Delete files if user is inactive
    fetch(`/delete/${localStorage.getItem("fileName")}`, { method: "DELETE" })
    console.log(`deleted file ${localStorage.getItem("fileName")}`)
  }

  const handleOnActive = event => {
    // Do nothing
  }

  const handleOnAction = event => {
    // Do nothing
  }

  const { getRemainingTime, getLastActiveTime } = useIdleTimer({
    timeout: 1000 * 60 * IDLE_MINUTES,
    onIdle: handleOnIdle,
    onActive: handleOnActive,
    onAction: handleOnAction,
    debounce: 500
  })

  return (
    <ThemeProvider theme={theme}>
      <CssBaseline />
      <Router>
        <div className="App">
          <Route exact path="/">
            <Upload />
          </Route>
          <Switch>
            <Route path="/menu">
              <Menu />
            </Route>
            <Route path="/visualize">
              <Visualize />
            </Route>
          </Switch>
        </div>
      </Router>
    </ThemeProvider>
  );
}

export default App;
