#import <Foundation/Foundation.h>

@interface BridgeClass : NSObject

-(NSDictionary<NSString*, NSNumber*>*)updateECG:(int)data;
-(NSDictionary<NSString*, NSNumber*>*)updatePPGRed:(double)Red IR:(double)IR;
@end
