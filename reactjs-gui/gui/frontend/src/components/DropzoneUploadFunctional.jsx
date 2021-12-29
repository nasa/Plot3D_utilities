// https://www.geeksforgeeks.org/file-uploading-in-react-js/
import axios from 'axios';
import React from 'react';
import { DropzoneDialog } from 'material-ui-dropzone'
import { Button, Input } from '@material-ui/core';
import UploadFileIcon from '@mui/icons-material/UploadFile';
import NavigateNextIcon from '@mui/icons-material/NavigateNext';
import { Link } from 'react-router-dom';


export default function DropzoneUploadFunctional() {
    // State
    const [open, setOpen] = React.useState(false);
    const [files, setFiles] = React.useState([]);
    const [filename, setFilename] = React.useState('');
    const [loading, setLoading] = React.useState(false);
    const [done, setDone] = React.useState(false);

    // Handlers
    const handleClose = () => {
        setOpen(false);
    };

    const handleSave = (files) => {
        setFiles(files);
        setOpen(false);
        setLoading(true);
        setDone(false);
        const formData = new FormData();
        formData.append('file', files[0]);
        axios.post('/upload', formData, {
            headers: {
                'Content-Type': 'multipart/form-data'
            }
        })
            .then(res => {
                localStorage.setItem('fileName', res.data.fileName);
                return res;
            })
            .then(res =>{
                setFilename(res.data.fileName);
                return res;
            })
            .then(res => {
                setLoading(false);
                setDone(true);
                return res;
            })
            .catch(err => {
                setLoading(false);
                setDone(false);
            });
    };

    const handleOpen = () => {
        setOpen(true);
    };

    const uploadFileButton = () => {
        if (loading) {
            return (
                <div>
                    Loading...
                </div>
            );
        } else if (done) {
            return (
                <div>
                    <Button
                        component={Link}
                        to="/visualize"
                        variant="contained"
                        color="primary"
                        startIcon={<NavigateNextIcon />}
                    >
                        Visualize your file
                    </Button>
                </div>
            );
        } else {
            return (
                <div>
                    <Button
                        variant="contained"
                        color="primary"
                        startIcon={<UploadFileIcon />}
                        onClick={handleOpen}
                    >
                        Upload File
                    </Button>
                    <DropzoneDialog
                        open={open}
                        onSave={handleSave}
                        onClose={handleClose}
                        acceptedFiles={[]}
                        showPreviews={true}
                        maxFileSize={1000000000}
                        filesLimit={1}
                        dropzoneText="Drag and drop an ASCII Plot3d file here (.xyz), or click to select a file to upload."
                        showPreviewsInDropzone={true}
                        showAlerts={true}
                        onChangeStatus={(file, status) => console.log(status)}
                        onSubmit={(files) => setFiles(files)}
                    />
                </div>
            );
        }
    }

    return(
        <>
            {uploadFileButton()}
        </>
    );
}