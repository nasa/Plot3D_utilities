import { useEffect, useState } from 'react';
import { Grid, Typography } from '@material-ui/core';

import DropzoneUploadFunctional from '../components/DropzoneUploadFunctional';

export default function Upload() {
  return (
    <Grid
      container
      spacing={0}
      direction="column"
      alignItems="center"
      justifyContent="center"
      style={{ minHeight: '100vh' }}
    >
      <Grid item xs={3}>
        <Typography variant="h4" gutterBottom>
          GlennHT GUI
        </Typography>
      </Grid>
      <Grid item xs={3}>
        <DropzoneUploadFunctional />
      </Grid>
    </Grid>
  );
}
