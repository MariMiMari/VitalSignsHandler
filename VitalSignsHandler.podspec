#
# Be sure to run `pod lib lint VitalSignsHandler.podspec' to ensure this is a
# valid spec before submitting.
#
# Any lines starting with a # are optional, but their use is encouraged
# To learn more about a Podspec see http://guides.cocoapods.org/syntax/podspec.html
#

Pod::Spec.new do |s|
  s.name             = 'VitalSignsHandler'
  s.version          = '0.0.1'
  s.summary          = 'VitalSignsHandler is ios framework. If you need to handler ECG or PPG data use it.'

# This description is used to generate tags and improve search results.
#   * Think: What does it do? Why did you write it? What is the focus?
#   * Try to keep it short, snappy and to the point.
#   * Write the description between the DESC delimiters below.
#   * Finally, don't worry about the indent, CocoaPods strips it!

  s.description      = <<-DESC
This framework was developer of Higher  School  of  Information  Technology  and  Information  Sys-tems of Kazan Federal University."
                       DESC

  s.homepage         = 'https://github.com/MariMiMari/VitalSignsHandler'
  # s.screenshots     = 'www.example.com/screenshots_1', 'www.example.com/screenshots_2'
  s.license          = { :type => 'MIT', :file => 'LICENSE' }
  s.author           = { 'Мария Тимофеева' => 'mi.maritim@gmail.com' }
  s.source           = { :git => 'https://github.com/MariMiMari/VitalSignsHandler.git', :tag => s.version.to_s }
  # s.social_media_url = 'https://twitter.com/<TWITTER_USERNAME>'

  s.ios.deployment_target = '8.0'

  s.libraries =  'c++'
  s.source_files = 'VitalSignsHandler/Classes/**/**' 
  # s.vendored_frameworks = 'VitalSignsHandler/VitalSignsHandler.framework' 
  # s.resource_bundles = {
  #   'VitalSignsHandler' => ['VitalSignsHandler/Assets/*.png']
  # }
  s.requires_arc = true
  s.xcconfig = {
     'CLANG_CXX_LANGUAGE_STANDARD' => 'c++11',
     'CLANG_CXX_LIBRARY' => 'libc++' ,
     # 'SWIFT_VERSION' => '4.0'
  }
  s.preserve_paths = 'Example/Pods/Target Support Files/VitalSignsHandler/VitalSignsHandler.modulemap'
  s.module_map = 'Example/Pods/Target Support Files/VitalSignsHandler/VitalSignsHandler.modulemap'

  s.public_header_files = 'VitalSignsHandler/Classes/**/*.h'
  # s.frameworks = 'UIKit', 'MapKit'
  # s.dependency 'AFNetworking', '~> 2.3'
end
