
Array.prototype.insert = function (z, item) {
	this.splice(z, 0, item);
};
Array.prototype.remove = function (z) {
	this.splice(z, 1);
};
Array.prototype.copy = function () {
	return this.slice();
};
String.prototype.nosuffix = function (x) {
	var z = this.lastIndexOf(x);
	if (z != -1)
		return this.substr(0, z);
	return this;
};

function uuid()
{
    var d = new Date().getTime();
    if(window.performance && typeof window.performance.now === "function"){
        d += performance.now(); //use high-precision timer if available
    }
    var uuid = 'xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx'.replace(/[xy]/g, function(c) {
        var r = (d + Math.random()*16)%16 | 0;
        d = Math.floor(d/16);
        return (c=='x' ? r : (r&0x3|0x8)).toString(16);
    });
    return uuid;
}

function api(func, arg, callback, aux = {}) 
{
	if(spmt.token != null)
		arg.token = spmt.token;

 	return $.ajax($.extend({url:spmt.apiurl+'/'+func, data:arg, dataType:'json', method: "POST", success:callback, async: true}, aux));
}

$(document).ready(
function ()
{
});
