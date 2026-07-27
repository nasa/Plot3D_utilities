// https://www.geeksforgeeks.org/file-uploading-in-react-js/
import axios from 'axios';
import React from 'react';
import { DropzoneDialog } from 'material-ui-dropzone'
import { Button, Input } from '@material-ui/core';
import UploadFileIcon from '@mui/icons-material/UploadFile';
import NavigateNextIcon from '@mui/icons-material/NavigateNext';
import { Link } from 'react-router-dom';

class DropzoneUpload extends React.Component {
    constructor(props) {
        super(props);
        this.state = {
            open: false,
            files: [],
            loading: false,
            doneUpload: false,
        };
    }


    handleClose() {
        this.setState({
            open: false
        });
    }
 
    handleSave(files) {
        const self = this;
        //Saving files to state for further use and closing Modal.
        this.setState({
            files: files,
            open: false,
            loading: true,
        });
        const formData = new FormData();

        formData.append(
            'file', files[0]
        );

        formData.append(
            "filename", files[0].name
        );

        axios.post("/upload", formData)
        .then(function (response) {
            console.log(response);
            self.setState({
                doneUpload: true,
                loading: false,
            });
        });

    }
 
    handleOpen() {
        this.setState({
            open: true,
        });
    }

    render() {
        const { match, location, history } = this.props

        if (!this.state.doneUpload) {
            return (
                <div>
                    <Button
                        variant="contained"
                        color="primary"
                        endIcon={<UploadFileIcon />}
                        onClick={this.handleOpen.bind(this)}
                    >
                    Upload File
                    </Button>
                    <DropzoneDialog
                        open={this.state.open}
                        onSave={this.handleSave.bind(this)}
                        acceptedFiles={[]}
                        filesLimit={1}
                        showPreviews={true}
                        maxFileSize={500000000}
                        onClose={this.handleClose.bind(this)}
                    />
                </div>
            );
        } else {
            if (this.state.loading){
                return (
                    <div>
                        Loading...
                    </div>
                );
            } else {
                return (
                    <div>
                        <Button
                            component={Link}
                            to="/menu"
                            variant="contained"
                            color="primary"
                            endIcon={<NavigateNextIcon />}
                        >
                        Next
                        </Button>
                    </div>
                );
            }
        }
    }
}

export default DropzoneUpload;