// https://blog.logrocket.com/react-suspense-for-data-fetching/
function WrapPromise(promise) {
    let status = "pending";
    let response

    const suspender = promise.then(
        (res) => {
            status = "success"
            response = res
        },
        (err) => {
            status = "error"
            response = err
        },
    )

    const read = () => {
        switch (status) {
            case "pending":
                throw suspender
            case "error":
                throw response
            default:
                return response
        }
    }

    return { read }
}

export default WrapPromise