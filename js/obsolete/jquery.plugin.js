(function ($) {
    // plugin 등록
    // Input: Class에 해당하는 function
    // Output: $의 함수의 chain기능을 할 수 있도록 
    //         $()의 jQuery object를 던져주는 것으로 기본으로 합니다.
    //        만일 값을 요구하는 것이면 
    //        jQuery object내 DOM Element가 한개이면 한개의 결과 값을 돌려주고
    //        여러개이면 array로 값을 돌려 줍니다.
	$.namespace = {};
    $.plugin = function (clsobj, parentObj) {
        var clsname = clsobj.name;
        if ($.fn[clsname]) return;
		if ( parentObj )
		{
			clsobj.__super = $.namespace[parentObj];
			clsobj.prototype = new clsobj.__super();
		}
		$.namespace[clsname] = clsobj;
        $.fn[clsname] = function (cmd) {
            var ret = this;
            var self = this;
            
            // $.clsname()에 입력 값이 없거나 {}나 new Object()의 값이면 생성을합니다.
            if ( !cmd || $.isPlainObject(cmd)) {
                $.each(this, function (idx, value) {
                    this.__instance = new clsobj(this, cmd);
                });
            }
            // 만일 첫번째 parameter가 string이면 함수명으로 인식하고 실행합니다.
            else if (typeof (cmd) == "string") {
                var params = Array.prototype.slice.call(arguments, 1);
                $.each(this, function (idx, value) {
                    if (this.__instance[cmd]) {
                        var v = this.__instance[cmd].apply(this.__instance, params);
                        if (v != undefined) {
                            if (self == ret) {
                                ret = [v];
                            }
                            else {
                                ret.push(v);
                            }
                        }
                    }
                });
            }
            if (ret == self) {
                return self;
            }
            return ret.length == 1 ? ret[0] : ret;
        }
    }
})(jQuery);
