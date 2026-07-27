import WrapPromise from "./WrapPromise";

const blockzeroURL = "/blockzero"
const blockoneURL = "/blockone"
const connectivity0URL = "/connectivity0"
const connectivity1URL = "/connectivity1"
const periodicityURL = "/periodicity"

function fetchBlockZeroURL() {
    const promise = fetch(blockzeroURL)
        .then((res) => res.json())
        .then((res) => res.data)

    return WrapPromise(promise)
}

function fetchBlockOneURL() {
    const promise = fetch(blockoneURL)
        .then((res) => res.json())
        .then((res) => res.data)

    return WrapPromise(promise)
}

function fetchConnectivity0URL() {
    const promise = fetch(connectivity0URL)
        .then((res) => res.json())
        .then((res) => res.data)

    return WrapPromise(promise)
}

function fetchConnectivity1URL() {
    const promise = fetch(connectivity1URL)
        .then((res) => res.json())
        .then((res) => res.data)

    return WrapPromise(promise)
}

function fetchPeriodicityURL() {
    const promise = fetch(periodicityURL)
        .then((res) => res.json())
        .then((res) => res.data)

    return WrapPromise(promise)
}

export { fetchBlockZeroURL, fetchBlockOneURL, fetchConnectivity0URL, fetchConnectivity1URL, fetchPeriodicityURL }