(function($){
    // UI관련 최상의 부모 Class
    function UIControl(element, properties) {
        // 상속시 아무런 parameter가 없기 때문에 이 경우 아무런 처리를 하지 않습니다.
        if (arguments.length == 0) return;
        var self = this;

        // 속성값들 저장
        this.element = $(element);
        if (this.properties) {
            this.properties = $.extend(this.properties, properties);
        }
        else {
            this.properties = properties || { enabled: true };
        }

        // 자식 생성
        var list = this.element.find('[data-child-name]');
        this.children = {};
        for (var i = 0 ; i < list.length; i++) {
            var child = $(list[i]);
            var name = child.data('childName');
            // 만일 plugin적용이 필요한가
            if (child.data('childClass')) {
                child[child.data('childClass')]();
            }
            this.children[name] = child;
        }

        // UI가 enabled가 아니면 아무런 처리를 하지 못하도록 하기 위해서 등록을 합니다.
        this.element.click(function (e) {
            if (!self.properties.enabled) {
                e.stopImmediatePropagation();
                e.stopPropagation();
                return false;
            }
        });
    }
    // plugin 등록
    $.plugin(UIControl);

    // 호출시 입력값이 없다면 현재 enabled값을 돌려 줍니다.
    // 부모를 disable하면 자식 모두 event를 받지 못하게 하지만,
    // 주의 할 점은 부모가 enable되면 모든 children도 enable됩니다.
    UIControl.prototype.enabled = function (value) {
        if (arguments.length > 0) {
            this.properties.enabled = value;
            for (var cname in this.children) {
                var child = this.children[cname];
                if (child[0].__instance) {
                    child.UIControl('enabled', value);
                }
            }
            if (!value) {
                this.element.addClass('disabled');
            }
            else {
                this.element.removeClass('disabled');
            }
        }
        return this.properties.enabled;
    }
})(jQuery);
