// https://www.geeksforgeeks.org/file-uploading-in-react-js/
import axios from 'axios';
import React from 'react';
import { Button, Input } from '@material-ui/core';

class FileUploadAxios extends React.Component {
    state = {
        // No file selected initially
        selectedFile: null,
        responseFile: null,
        
    };

    // On file select
    onFileChange = event => {
        // Update state
        this.setState( {selectedFile: event.target.files[0] });
    };

    // On file upload
    onFileUpload = () => {
        // Create FormData object
        const formData = new FormData();

        // Update FormData object
        formData.append(
            "file", this.state.selectedFile,
        );

        formData.append(
            "filename", this.state.selectedFile.name
        );

        // Details of the uploaded file
        console.log(this.state.selectedFile);

        // Request made to backend api, send FormData object
        axios.post("/upload", formData)
            .then(function (response) {
                console.log(response);

        });
    };

    // File content displayed after file upload complete
    fileData = () => {
        if (this.state.selectedFile) {
            return (
                <div>
                    <h4>Click Confirm Upload if the following is correct:</h4>
                    <p>File Name: {this.state.selectedFile.name}</p>
                    <p>File MIME Type: {this.state.selectedFile.type}</p>
                    <p>File Last Modified: {new Date(this.state.selectedFile.lastModified).toISOString()}</p>
                </div>
            );
        }
        else {
            return (
                <div>
                    <br />
                    <h4>Upload a file by clicking Choose File</h4>
                </div>
            );
        }
    };

    render() {
        return (
            <div>
                <h1>File Upload</h1>
                <div>
                    <Input type="file" onChange={this.onFileChange} />
                    <Button variant="contained" color="primary" onClick={this.onFileUpload}>
                        Confirm Upload
                    </Button>
                </div>
                {this.fileData()}
            </div>
        );
    }
}

export default FileUploadAxios;