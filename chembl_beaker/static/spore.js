// SPoRE - jQuery port
// niko@rtgi


Raphael.st.draggable = function() {
  var me = this,
      lx = 0,
      ly = 0,
      ox = 0,
      oy = 0,
      moveFnc = function(dx, dy) {
          lx = dx + ox;
          ly = dy + oy;
          me.transform('t' + lx + ',' + ly);
      },
      startFnc = function() {},
      endFnc = function() {
          ox = lx;
          oy = ly;
      };

  this.drag(moveFnc, startFnc, endFnc);
};

//----------------------------------------------------------------------------------------------------------------------

var spore = {

//----------------------------------------------------------------------------------------------------------------------

    fileToUpload : null,

    handleFileSelect: function (evt) {
        var target = $(evt.currentTarget);
        evt.dataTransfer = evt.originalEvent.dataTransfer;
        evt.stopPropagation();
        evt.preventDefault();

        var files = evt.dataTransfer? evt.dataTransfer.files : evt.target.files; // FileList object.
        var file = files[0];
        var name = file.name;

        if (!name.match('(.*)\\.(jpg|gif|png|bmp|tiff|pdf|mol|sdf|mrv|txt|xml|json)'))
            return alert('Unrecognised file extension.');

        if(files.length > 1)
            return alert('More than one file added', 'Sorry, you can currently add only one file at time...');

        if(file.size > 2000000)
            return alert('File too large', 'Files larger than 2MB will not be accepted now.');

        var descriptionBox = $('.dropDesc', target);
        descriptionBox.text(name);
        var reader = new FileReader();
        reader.onload = function(e){
            spore.fileToUpload = e.target.result;
        };

        reader.onerror = function(error){
            var desc  = error.target.error.name ? error.target.error.name : '';
            alert("Error reading file: " + desc);
        };

        if (name.match('(.*)\\.(jpg|gif|png|bmp|tiff|pdf)'))
            reader.readAsDataURL(file);
        else if(name.match('(.*)\\.(mol|sdf|mrv|txt|xml|json)'))
            reader.readAsText(file);

    },

//----------------------------------------------------------------------------------------------------------------------

    handleDragOver: function (evt) {
        evt.dataTransfer = evt.originalEvent.dataTransfer;
        evt.stopPropagation();
        evt.preventDefault();
        evt.dataTransfer.dropEffect = 'copy'; // Explicitly show this is a copy.
    },

//----------------------------------------------------------------------------------------------------------------------

    handleClipboard: function (evt) {
        var response = $('.rawResponse', $(evt.target).closest('.restapiholderresponse')).text();
        ClipboardCopyInit(this, response);
    },

//----------------------------------------------------------------------------------------------------------------------

    log: function () {},

//----------------------------------------------------------------------------------------------------------------------

    debug: function (state) {
        if (state) {
            var pfx = '[spore] ';
            if (window.console && window.console.debug) //firebug
            spore.log = function (m) {
                window.console.debug(pfx + m);
            };
            else if (window.console && window.console.log) //safari
            spore.log = function (m) {
                window.console.log(pfx + m);
            };
            else if (console && console.log) //chrome
            spore.log = function (m) {
                console.log(pfx + m);
            };
            else if (opera && opera.postError) //opera
            spore.log = function (m) {
                opera.postError(pfx + m);
            };
        } else spore.log = function () {};
    },

//----------------------------------------------------------------------------------------------------------------------

    create: function (spec, callback, over,destination) {
        if (typeof (spec) == 'string') {
            if (/^\s*{/.test(spec)) spec = $.parseJSON(spec);
            else {
                //todo: spec cache
                if (!callback) throw "your must provide a callback with a spec over http";
                $.get(spec, function (data) {
                    spore.create(data, callback, over,destination)
                });
                return true;
            }
        } else if (typeof (spec) != 'object') throw "wrong spore creation parameters";

        if (over) spec = $.extend(spec, over);
        var a = new spore.api(spec,'',destination);
        if (callback) callback(a);
        return a;
    },

//----------------------------------------------------------------------------------------------------------------------

    api: function (spec, widgets,destination) {
        this.spec = $.extend({
            api_format: 'json',
            authentication: false,
            name: 'unknown',
            version: '?.?'
        }, spec);
        this.widgets = [];
        this.destination = destination; //This the destination we're going to build the widgets in (IODOC stylee)

        spore.log('spore init api ' + this.spec.name + ' v' + this.spec.version);

        //create methods
        var self = this;
        var wrapfn = function (fn) {
            return function (params, onsuccess, onerror) {
                self._call(fn, params, onsuccess, onerror);
            };
        };

        $(document).on('dragover', '.drop_zone', spore.handleDragOver);
        $(document).on('drop', '.drop_zone', spore.handleFileSelect);
        $(document).on('change', '.fileUploader', spore.handleFileSelect);
        $(document).on('mouseover', '.clipboardBtn', spore.handleClipboard);

        spore.settitle(this.destination,this.spec);


        for (var fn in this.spec.methods) {
            if (!/^\//.test(this.spec.methods[fn].path)) this.spec.methods[fn].path = '/' + this.spec.methods[fn].path;
            this[fn] = wrapfn(fn);
            spore.log(this.spec.name + '.' + fn + ' created');

            spore.createwidget(this.destination,fn);

        }
    },

//----------------------------------------------------------------------------------------------------------------------

    settitle: function(destination,spec) {
        this.destination = destination;
        this.spec = spec;

        $('#'+this.destination).append('<div class="span12" id="mainholder"></div>');

        var holder = $('#mainholder');
        holder.append('<div id="title"><h1>'+this.spec.name + ' Explorer</h1><p>Version: ' +
                                                                                this.spec.version+'</p></div>');
        holder.append('<div>Your (optional) API Key: <input type="text" name="chembl_api_key" id="chembl_api_key"></div>');
        holder.append('<div class="accordion" id="restaccordian"></div>');
    },

//----------------------------------------------------------------------------------------------------------------------

    createwidget: function(destination,fn) {

        // Use jQuery and Bootstrap to create a widget to test the rest endpoint.

        var params={};
        var apifunction = this.spec.methods[fn];
        var method = apifunction.method;
        var path = apifunction.path;
        var description = apifunction.description;
        var required_params = apifunction.required_params;
        var formats = apifunction.formats;

        params[required_params]='';

        var splitter = ':'+required_params;
        var testtailparts = path.split(splitter);
        var testtail=testtailparts[1];
        var testhead=testtailparts[0];
        var btnclass='';
        var iconclass='';
        var imageUpload = (path.indexOf("image2") != -1);

        if(method=='GET')
            {
                btnclass='btn btn-success ';
                iconclass='icon-circle-arrow-down icon-black';
            }

        else if(method=='PUT' || method=='POST')
            {
                btnclass='btn btn-warning';
                iconclass='icon-circle-arrow-up icon-black';

            }

        else if(method=='HEAD')
            {
                //Ommiting HEAD requests
                return;
                //btnclass='btn btn-info';
                //iconclass='icon-ok-sign icon-black';
            }

        else if(method=='DELETE')
            {
                btnclass='btn btn-danger';
                iconclass='icon-remove-sign icon-black';
            }

        else if(method=='OPTIONS')
            {
                return;
            }

        var pathcopy = path;

        if(required_params) {

                $.each(required_params,function(key,value){

                    var re = new RegExp(':'+value,"g");
                    if(value == "CTAB")
                        newpart = '<textarea id="'+fn+'_input_'+value+'" class="input-small" type="text" placeholder="'+value+'"/>';
                    else if (imageUpload && value=='IMAGE')
                        newpart = '<div class="drop_zone well" id="'+fn+'_file_upload_'+ method + '"><div class="dropDesc">Drop image here</div></div>';
                    else
                        newpart = '<input id="'+fn+'_input_'+value+'" class="input-small" type="text" placeholder="'+value+'"/>';

                    pathcopy = pathcopy.replace(re, newpart);

                });
        }

        var display_path = path;
        if(display_path[0] == '/')
            display_path = display_path.substring(1);
        $('#restaccordian').append('<div class="accordion-group"><div class="accordion-heading"><a class="accordion-toggle" data-toggle="collapse" data-parent="#restaccordian" href="#collapse_'+fn+'"><button style="display:inline-block" class="'+btnclass+'">'+method+'</button><p style="display:inline-block;font-size:large;padding-left:10px">'+display_path+'</p></a></div><div id="collapse_'+fn+'" class="accordion-body collapse in"><div class="accordion-inner '+fn+'"></div></div></div>');

        var holder_widget = $('<div class="restapiholder" id="'+fn+'_holder">');
        var method_widget = $('<div class="restapiholdermethod" id="'+fn+'_method"></div>');
        $("."+fn).append(holder_widget);

        holder_widget.append(method_widget);
        holder_widget.append('<div class="restapiholderresponse" id="'+fn+'_response"></div>');
        method_widget.append('<p><h4>Description</h4><div class="well" id="'+fn+'_method_description"></div></p>');

        var base_url = (this.spec.base_url != undefined) ? this.spec.base_url : window.location.href.substr(0, window.location.href.lastIndexOf('/')+1);
        description = description.replace(/\$\{BEAKER_ROOT_URL\}/g, base_url);
        if(navigator.platform.match(/(Mac|iPhone|iPod|iPad)/i)){
            description = description.replace(/-w 0/g, '-b 0');
        }

        $('#'+fn+'_method_description').html(markdown.toHTML(description));

        if(required_params != undefined){
            method_widget.append('<p><h4>Requires</h4><div class="well">'+required_params+'</div></p>');
        } else {
            method_widget.append('<p><h4>Requires</h4><div class="well"><italic>No parameters required</italic></div></p>');
        }


        method_widget.append('<p><h4>Formats</h4><div class="well">'+formats+'</div></p>');
        if(required_params != undefined)
        {
            if (method=='POST') {
                if(!imageUpload)
                    method_widget.append('<textarea id="'+fn+'_input_post" rows="3" placeholder="POST body"></textarea>');
                else
                    method_widget.append('<div class="drop_zone well" id="'+fn+'_file_upload_'+ method + '"><div class="dropDesc">Drop image here</div></div>');
                method_widget.append('...or upload a file: <div><form style="margin-top: 1em;"><input type="file" class="fileUploader" name="file[]" multiple /></form></div>');
                required_params.push('post');
            }

            method_widget.append('<p>Enter a value for <italic>'+required_params+'</italic> and click '+method+' to test the service!</p>');
            method_widget.append('<p><button id="'+fn+'_trigger" class="'+btnclass+'"><i class="'+iconclass+'"></i> '+method+'</button><strong> '+pathcopy);
          } else {

            if (method=='POST') {
                if(!imageUpload)
                    method_widget.append('<textarea id="'+fn+'_input_post" rows="3" placeholder="POST body"></textarea>');
                else
                    method_widget.append('<div class="drop_zone well" id="'+fn+'_file_upload_'+ method + '"><div class="dropDesc">Drop image here</div></div>');
                method_widget.append('...or upload a file: <div><form style="margin-top: 1em;"><input type="file" class="fileUploader" name="file[]" multiple /></form></div>');
                var required_params = [];
                required_params.push('post');
            }

            method_widget.append('<p>Click '+method+' to test the service!</p>');
            method_widget.append('<p><button id="'+fn+'_trigger" class="'+btnclass+'"><i class="'+iconclass+'"></i> '+method+'</button>');

        }
    $('#collapse_'+fn).collapse("hide");


    $(document).ajaxStart(function () {

        $('#'+fn+'_response').html('<div id="'+fn+'_progress" class="restloadingDiv"></div>');
             }).ajaxStop(function () {

                    $('#'+fn+'_progress').remove();
                      $('#'+fn+'_highlight').each(function(i, e) {
                                hljs.highlightBlock(e, null, false)
                         });
        });

    var triggerflags;
    var run;

    $('#'+fn+'_trigger').click({required_params:required_params},function(event){

       triggerflags = [];
       run = true;

       if(spore.fileToUpload && method == 'POST') {
           params['post']=spore.fileToUpload;
       }

       else if(required_params) {
        $.each(required_params,function(key,value){
            var input = $('#'+fn+'_input_'+value);
            if($.inArray(value, ['CTAB', 'SMILES', 'INCHI', 'SMARTS', 'TEMPLATE']) != -1 && method == 'GET'){
                params[value]= $.base64.urlsafe_encode(input.val());
            }
            else if(value == 'IMAGE'){
                params[value] = spore.fileToUpload.replace(/\+/g, '-').replace(/\//g, '_');
            }
            else{
                params[value]= input.val();
            }

            if(params[value] !=undefined){
                if(params[value].length==0){

                     input.popover({html:false,title:'Error',content:'No value entered!',delay: { show: 0, hide: 100 },trigger:'focus'});
                     input.focus();

                     triggerflags.push(false);

                } else {

                     triggerflags.push(true);

                     input.popover('disable');
                }
             }
        });

        $.each(triggerflags,function(key,value){
            // If we have a false flag, i.e. no argument => return
            if(!value){
                    run = false;
                }
            });

       if(!run)
           return;
       }

       spore.fileToUpload = null;

    spore.log("RP: "+required_params);

        if (params['post']!=undefined && required_params.length == 1) {
            api.spec.methods[fn].required_params=['post'];
            spore.log(api.spec.methods[fn]);
        }

        spore.log('params: '+JSON.stringify(params)+' '+method);

        api._call(fn,params,

                 function(data, req, url){
                    // Success (jQuery is misbehaving, switching to native JS)
                    var content_type = req.getResponseHeader('Content-Type');
                    document.getElementById(fn+'_response').classList.remove('resterror');
                    document.getElementById(fn+'_response').classList.add('restsuccess');
                    var raw_data = '';
                    if($.type(data) === "string")
                        raw_data = data;
                    else if ($.isXMLDoc(data)){
                        raw_data = (new XMLSerializer()).serializeToString(data);
                    }
                    else
                        raw_data = JSON.stringify(data);
                    $('#'+fn+'_response').html('<h4>Request URI</h4><div id="'+fn+'_requestURI" class="mhsy well">'+url+'</div><h4>Response Code</h4><div class="well">'+req.status+'</div><h4>Response</h4><div class="well">'+req.statusText+'</div><h4>Response Body <button class="clipboardBtn">Copy</button><span class="clipboardStatus" /></h4><div id="'+fn+'_responseJSON" class="responseBody"></div>');
                    var raw_response = $('<div class="rawResponse" style="display:none;"></div>');
                    raw_response.text(raw_data);
                    $('#'+fn+'_response').append(raw_response);
                    console.log('FN = ' + fn);
                    if(fn.indexOf("23D") != -1){
                            add_molecule_canvas('', '', raw_data, $('#'+fn+'_responseJSON'));
                    }
                    else if (content_type.indexOf('image/png') != -1){
                        $('#'+fn+'_responseJSON').html('<img src="data:image/png;base64,' + data + '" />');
                    }
                    else if(content_type.indexOf('application/json') != -1){
                        var cell_size = parseInt(params['SIZE']) || 200;
                        var canvas = $('<div/>').width(cell_size).height(cell_size);
                        $('#'+fn+'_responseJSON').append(canvas);
                        var paper = Raphael(canvas[0], cell_size, cell_size);
                        var compound = paper.add(data);
                        compound.draggable();
                        var curentSize = cell_size;
                        canvas.bind('mousewheel DOMMouseScroll', function(event) {
                                            event.preventDefault();
                                            var delta = event.originalEvent.wheelDelta || -event.originalEvent.detail;
                                            curentSize *= (1.0 + delta / 100.0);
                                            var x = (paper.width / 2) - (curentSize /2);
                                            var y = (paper.height / 2) - (curentSize /2);
                                            paper.setViewBox(x, y, curentSize, curentSize, false);

                        });
                    }
                    else if (content_type.indexOf('image/svg+xml') != -1){
                        var importedSVGRootElement = document.importNode(data.documentElement,true);
                        $('#'+fn+'_responseJSON').html(importedSVGRootElement);
                    }
                    else if (data[0] == '{' || data[0] == '['){
                        try{
                           $('#'+fn+'_responseJSON').html('<pre id="'+fn+'_pre"><code id="'+fn+'_highlight">'+vkbeautify.json(data)+'</code></pre>');
                        }
                        catch(err){
                           $('#'+fn+'_responseJSON').html('<pre id="'+fn+'_pre"><code id="'+fn+'_highlight">' + data + '</code></pre>');
                        }
                    }

                    else{
                        $('#'+fn+'_responseJSON').html('<pre id="'+fn+'_pre"><code id="'+fn+'_highlight">'+ data +'</code></pre>');
                    }

                    $('#'+fn+'_highlight').addClass('mhsy');
                    $('#'+fn+'_requestURI').addClass('mhsy');


                                        },
                 function(req,status,err,url){
                    // Error
                    try{
                        errorobj=JSON.parse(req.responseText);
                    }
                    catch(err){
                        spore.log(JSON.stringify(req));
                        // Split out the request for 200 (seperate the error/sucess methods)

                        errorobj=err;
                    }

                    uri = errorobj.request_uri;
                    error = errorobj.error;
                    document.getElementById(fn+'_response').classList.remove('restsuccess');
                    document.getElementById(fn+'_response').classList.add('resterror');
                    document.getElementById(fn+'_response').innerHTML = '<h4>Request URI</h4><div class="mhsy well">'+uri+'</div><h4>Response Code</h4><div class="well">'+req.status+'</div><h4>Response</h4><div class="well">'+req.statusText+'</div><h4>Error Body</h4><div id="'+fn+'_responseJSON"></div>';
                    $('#'+fn+'_responseJSON').html('<pre id="'+fn+'_pre"><code id="'+fn+'_highlight">'+JSON.stringify(error,undefined,2)+'</code></pre>');

                    $('#'+fn+'_highlight').addClass('mhsy');
                    $('#'+fn+'_requestURI').addClass('mhsy');
                 }
            );
    });

    },
    widgets: {}
};

//----------------------------------------------------------------------------------------------------------------------

spore.api.prototype._call = function (fn, params, onsuccess, onerror) {
    var method = this.spec.methods[fn],
        url = (this.spec.base_url != undefined) ? this.spec.base_url + method.path : method.path.substring(1),
        prm = $.extend({}, params),
        plog = this.spec.name + '.' + fn + ' ';
    dta = {};
    spore.log(plog + 'init(' + (params == undefined ? '' : (JSON ? JSON.stringify(params) : '...')) + ')');
    //manage params
    if (method.required_params) for (var n = 0; n < method.required_params.length; ++n) {
        var p = method.required_params[n];
        if (prm[p] == undefined) throw 'param ' + p + ' required';
        var r = new RegExp("(\:" + p + ")(\/|$)");
        if (r.test(url)) url = url.replace(r, prm[p] + '$2');
        else dta[p] = prm[p];
    }
    if (method.optional_params) for (var n = 0; n < method.optional_params.length; ++n) {
        var p = method.optional_params[n],
            d = (prm[p] != undefined),
            r = new RegExp("(\:" + p + ")(\/|$)");
        if (r.test(url)) url = url.replace(r, d ? prm[p] + '$2' : '$2');
        else if (d) dta[p] = prm[p];
    }
    //call widgets
    for (var n = 0; n < this.widgets.length; n++) {
        var w = this.widgets[n],
            wlog = plog + 'widget ' + w.name + ' ';
        try {
            var o = w.call(req, dta);
            if (o === false) {
                spore.log(wlog + 'returned false, call stopped!');
                return false;
            } else if (typeof (o) == 'function') {
                spore.log(wlog + 'returned a callback');
                w.callback = o;
            } else {
                spore.log(wlog + 'returned ' + o);
            }
        } catch (e) {
            spore.log(wlog + 'error! => ' + e);
        }
    }
    //webservice's call
    spore.log(plog + 'call ' + url);
    var self = this;
    dta['format'] = 'json';

    if (method.method=='POST') {
        dta=dta['post'];
    }

    var headers = {'X-Requested-With': 'XMLHttpRequest'};
    var api_key = $.trim($('#chembl_api_key').val());
    if(api_key)
        headers['X-ChEMBL-APIKey'] = api_key;

    $.ajax({
        url: url,
        type: method.method,
        data: dta,
		crossDomain: false,
		headers: headers,
        //dataType: 'text',
        //contentType: 'json', //this.spec.formats[0]
        success: function (data, status, req) {
            spore.log(plog + 'return ' + (data ? "datas" : "nothing"));
            //call widgets backwards
            for (var n = self.widgets.length - 1; n >= 0; --n) {
                var w = self.widgets[n],
                    wlog = plog + 'widget ' + w.name + ' ';
                if (w.callback) {
                    try {
                        var r = w.callback(data, req);
                        spore.log(wlog + 'callback returned ' + r);
                    } catch (e) {
                        spore.log(wlog + 'callback error! => ' + e);
                    }
                    delete w.callback;
                }
            }
            if (onsuccess) onsuccess(data,req,url);
        },
        error: function (req, status, err, url) {
            spore.log(plog + 'server returned '+url+' '+ req.status + ':' + req.statusText);
            if (onerror) onerror(req,status,err,url);
        }
    });
    return true;
};

//----------------------------------------------------------------------------------------------------------------------

spore.api.prototype.enable = function (widgetname, params) {
    if (!spore.widgets[widgetname]) throw 'widget ' + widgetname + ' not found/loaded';
    this.widgets.push({
        name: widgetname,
        obj: new spore.widgets[widgetname](params)
    });
    return this;
};

//----------------------------------------------------------------------------------------------------------------------